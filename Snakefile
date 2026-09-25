include: "workflows/test.smk"

rule default:
    default_target: True
    input:
        rules.all.input,
