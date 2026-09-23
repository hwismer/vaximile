process PREPARE_FASTA {

    // Unzips a fasta.gz or fa.gz or otherwise renames to a standard name

	label 'process_low'

    tag "Preprocessing $fasta"

    input:
    	path fasta

    output:
    	path "*_prc.fa", emit: fasta

    script:

    	def is_gz = fasta.name.endsWith('.gz')
        def prefix = task.ext.prefix ?: fasta.name.replaceFirst(/\.(fasta|fa)(\.gz)?$/, '')

    	"""
    	set -euo pipefail

    	if ${is_gz}; then
        	gunzip -c ${fasta} > ${prefix}_prc.fa
    	else
        	cp ${fasta} ${prefix}_prc.fa
    	fi
    	"""
    stub:

        def prefix = task.ext.prefix ?: fasta.name.replaceFirst(/\.(fasta|fa)(\.gz)?$/, '')

        """
        touch ${prefix}_prc.fa
        """

}
