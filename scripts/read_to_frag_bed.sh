#########1#########2#########3#########4#########5#########6#########7#########8

# Snakemake variables
input_bam="$1"
params_fasta="$2"
output_frag_bed="$3"

# Function
bam_to_frag(){
    # Make bedpe
    bedtools bamtobed -bedpe -i $1 |
        # Filter any potential non-standard alignments
        awk '$1==$4 {print $0}' | awk '$2 < $6 {print $0}' |
        # Create full-fragment bed file
        awk -v OFS='\t' '{print $1,$2,$6}' |
        # Annotate with GC content and fragment length
        bedtools nuc -fi $2 -bed stdin |
        # Convert back to standard bed with additional columns
        awk -v OFS='\t' '{print $1,$2,$3,$5,$12}' |
        sed '1d' > $3
    }

# Run command
bam_to_frag $input_bam \
            $params_fasta \
            $output_frag_bed
