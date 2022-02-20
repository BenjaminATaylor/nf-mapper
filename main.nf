#!/usr/bin/env nextflow

// Configurable variables
params.publish_dir = '../../../output'
samplelist = params.samplelist   
rawdatadir =  "realpath $params.rawdatadir".execute().text.trim()
genome = params.genome

// Validate inputs
if ( samplelist.isEmpty () ) {
    exit 1, 'Please specify --samplelist with a list of expected sample IDS'
} else if ( rawdatadir.isEmpty () ) {
    exit 1, 'Please specify --rawdatadir with the directory containing raw data'
}

//pull the list of samples and assign to an array
File samplels = new File(samplelist)
def samples = samplels.readLines()
//list files in the raw data directory
File rawls = new File(rawdatadir)
rawfiles = rawls.listFiles()
//check that each sample ID corresponds to exactly two input fastas
def countFastas(ID) {
    def pattern = ID + "_[1-2].fq.gz"
    def matcher = rawfiles =~ pattern
    return matcher.count
}
if ( !samples.every(it -> countFastas(it) == 2) ) {
    exit 1, 'Every supplied sample ID should correspond to exactly two fastq files.'
}

samples_in_trimmomatic = Channel.fromList(samples)
rawdatadirchan = Channel.value( rawdatadir )
genomechan = Channel.value( genome )

process trimmomatic {

    publishDir '${params.publish_dir}/trimmomatic', mode: 'copy'
    module 'bioinfo:trimmomatic/0.39'
    clusterOptions '--ntasks 16 --time 30:00 -A bharpur'

    input:
    val sampleID from samples_in_trimmomatic
    val rawpath from rawdatadirchan

    output:
    tuple val(outid), val(sampleID) into ch_out_trimmomatic

    script:

    fq_1 = rawpath + '/' + sampleID + '_1.fq.gz'
    fq_2 = rawpath + '/' + sampleID + '_2.fq.gz'
    outdir = params.publish_dir + '/trimmomatic/' + sampleID
    outid = outdir + '/' + sampleID + ".fq.gz"

    //TODO: re-parameterize the trimmomatic command for optimal trimming
    """
    mkdir -p $outdir

    trimmomatic \
    PE \
    $fq_1 \
    $fq_2 \
    -baseout $outid \
    -threads 16 \
    LEADING:3 TRAILING:3 MAXINFO:36:0.7 MINLEN:36
    """
}

process nextgenmap{

    publishDir '{params.publish_dir}/nextgenmap', mode: 'copy'
    module 'bioinfo:samtools/1.12'
    memory '10 GB'
    clusterOptions '--ntasks 16 --time 1:00:00 -A bharpur'

    input:
    tuple val(trimpath), val(sampleID) from ch_out_trimmomatic
    path genome from genomechan

    output: 
    tuple val(outfile), val(sampleID) into ch_out_nextgenmap

    script:
    trimfile1 = trimpath.replaceAll(".fq", "_1P.fq")
    trimfile2 = trimpath.replaceAll(".fq", "_2P.fq")

    outdir = params.publish_dir + '/nextgenmap/' + sampleID
    outfile = outdir + "/" + sampleID + "_ngm_sorted.bam"

    """
    rm -rf $outdir; mkdir -p $outdir

    ngm -r $genome \
    --qry1 $trimfile1 --qry2 $trimfile2 \
    -t 16 | 
    samtools view -bu | 
    samtools sort -o $outfile
    """

}

process qualimap{

    publishDir '{params.publish_dir}/qualimap', mode: 'copy'
    module 'bioinfo:qualimap/2.2.1'
    clusterOptions '--ntasks 16 --time 6:00:00 -A bharpur'

    input:
    tuple val(ngmout), val(sampleID) from ch_out_nextgenmap

    output:
    val outdir into ch_out_qualimap

    script: 
    outdir = params.publish_dir + '/qualimap/' + sampleID 
    """
    mkdir -p outdir
    unset DISPLAY

    qualimap bamqc -bam $ngmout --java-mem-size=2G -outdir $outdir
    """
}

process qualimap_prep{

    module 'r'
    clusterOptions '--time 10:00 -A bharpur'

    input:
    val qualdirs from ch_out_qualimap.collect()

    output:
    file 'outframe.tsv' into ch_out_qualprep

    script:
    qualdirString  = "\"" + qualdirs + "\""

    """
    #!/apps/spack/bell/apps/r/4.1.2-gcc-9.3.0-rw7vp7m/bin/Rscript

    qualdirString = gsub("]","",gsub("[","",$qualdirString,fixed=TRUE), fixed = T) 
    qualdirString = unlist(strsplit(qualdirString, split = ","))

    sampIDstring = gsub(".*qualimap/","", qualdirString)

    write.table(data.frame(cbind(sampIDstring, qualdirString)),file="outframe.tsv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
    """
}

process qualimap_collate{
    
    publishDir '{params.publish_dir}/qualimap_multiqc', mode: 'copy'
    module 'bioinfo:qualimap/2.2.1'

    input:
    file qualframe from ch_out_qualprep

    script:
    outdir = params.publish_dir + "/qualimap_multiqc/"

    """
    mkdir -p $outdir
    qualimap multi-bamqc -d $qualframe -outdir $outdir
    """
}

