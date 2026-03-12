#!/usr/bin/env python3
"""Post-process org HTML export to inline CSS and images as base64."""
import base64
import re
import sys
from pathlib import Path

def inline_html(html_path, project_root=None):
    html_path = Path(html_path)
    html_dir = html_path.parent
    project_root = Path(project_root) if project_root else html_dir
    html = html_path.read_text()

    # Inline CSS <link> tags
    def replace_css(m):
        href = m.group(1)
        css_path = (html_dir / href).resolve()
        if css_path.exists():
            css = css_path.read_text()
            return f'<style type="text/css">{css}</style>'
        return m.group(0)

    html = re.sub(
        r'<link\s+rel="stylesheet"\s+type="text/css"\s+href="([^"]+)"\s*/?>',
        replace_css, html)

    # Inline <img> tags
    def replace_img(m):
        full_tag = m.group(0)
        src = m.group(1)
        if src.startswith(('http://', 'https://', 'data:')):
            return full_tag
        img_path = (html_dir / src).resolve()
        if not img_path.exists():
            img_path = (project_root / src).resolve()
        if img_path.exists():
            mime = 'image/png' if img_path.suffix == '.png' else 'image/jpeg'
            b64 = base64.b64encode(img_path.read_bytes()).decode()
            return full_tag.replace(src, f'data:{mime};base64,{b64}')
        return full_tag

    html = re.sub(r'<img\s[^>]*src="([^"]+)"[^>]*/?>',  replace_img, html)

    html_path.write_text(html)
    print(f"Inlined: {html_path} ({html_path.stat().st_size // 1024} KB)")

if __name__ == '__main__':
    root = sys.argv[2] if len(sys.argv) > 2 else None
    inline_html(sys.argv[1], root)
