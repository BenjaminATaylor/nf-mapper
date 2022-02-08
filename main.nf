#!/usr/bin/env nextflow

// Configurable variables
def samplelist = '/depot/bharpur/data/popgenomes/USDAHornets/trimmomatic/hornet_IDs.txt'
def masterdir = '/depot/bharpur/data/popgenomes/USDAHornets/raw_data'
def genome = '/depot/bharpur/data/ref_genomes/VMAN/GCF_014083535.2_V.mandarinia_Nanaimo_p1.0_genomic.fna'

// Validate inputs
if ( samplelist.isEmpty () ) {
    exit 1, 'Please specify --samplelist with a list of expected sample IDS'
} else if ( masterdir.isEmpty () ) {
    exit 1, 'Please specify --masterdir with the directory containing raw data directories'
}

//for testing purposes
println 'Hello World'

//pull the list of samples and assign to an array
File samplels = new File(samplelist)
def samples = samplels.readLines()
//println samples

//for testing purposes, reassign samples as a single sample
samples = [samples.get(0), samples.get(1)]
println samples

//check that each sample has a corresponding directory in master
masterdirs = "ls $masterdir".execute().text
def mastercontains = { it -> masterdirs.contains(it) }
if ( !samples.every(mastercontains) ) {
    //println "Every sample in samplelist should have a corresponding folder in the master directory."
    exit 1, 'Every sample in samplelist should have a corresponding folder in the master directory.'
}

//check that each sample directory contains exactly two fastq files
def countfastas = { it -> new File("$masterdir/$it").listFiles() .findAll { it.name =~ /(?i)fq|fastq/ } .size() }
if ( !samples.every(it -> countfastas(it) == 2) ) {
    //println "Every sample directory should contain exactly two fastq files."
    exit 1, 'Every sample directory should contain exactly two fastq files.'
}

samples_in_trimmomatic = Channel.fromList(samples)
masterdirchan = Channel.value( masterdir )
genomechan = Channel.value( genome )

process trimmomatic {
    
    module 'bioinfo:trimmomatic/0.39'

    input:
    val sampleID from samples_in_trimmomatic
    path masterdir from masterdirchan

    output:
    val outid into ch_out_trimmomatic

    script:

    fq_1 = masterdir + '/' + sampleID + '/*_1.fq.gz'
    fq_2 = masterdir + '/' + sampleID + '/*_2.fq.gz'
    //fq_1_paired = '/depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/trimmomatic/' + sampleID + '/' + sampleID + '_1_trim_paired.fq.gz'
    //fq_1_unpaired = '/depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/trimmomatic/' + sampleID + '/' + sampleID + '_1_trim_unpaired.fq.gz'
    //fq_2_paired = '/depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/trimmomatic/' + sampleID + '/' + sampleID + '_2_trim_paired.fq.gz'
    //fq_2_unpaired = '/depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/trimmomatic/' + sampleID + '/' + sampleID + '_2_trim_unpaired.fq.gz'
    outdir = '/depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/trimmomatic/' + sampleID
    outid = outdir + '/' + sampleID + ".fq.gz"

    println outid

    //TODO: re-parameterize the trimmomatic command for optimal trimming
    """
    mkdir -p $outdir

    trimmomatic \
    PE \
    $fq_1 \
    $fq_2 \
    -baseout $outid \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

process nextgenmap{

    module 'bioinfo:samtools/1.12'
    memory '10 GB'
    clusterOptions '--ntasks 16 --time 30:00 -A bharpur'

    input:
    val trimpath from ch_out_trimmomatic
    path masterdir from masterdirchan
    path genome from genomechan

    output: 
    val outfile into ch_out_nextgenmap

    script:
    def sampID = (trimpath =~ /AGH[^\/]*/)[0]

    trimfile1 = trimpath.replaceAll(".fq", "_1P.fq")
    trimfile2 = trimpath.replaceAll(".fq", "_2P.fq")

    outdir = '/depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/nextgenmap/' + sampID
    outfile = outdir + "/" + sampID + "_ngm_sorted.bam"

    """
    echo $outfile
    echo $trimfile1
    echo TESTDONE

    mkdir -p $outdir

    ngm -r $genome \
    --qry1 $trimfile1 --qry2 $trimfile2 \
    -t 16 | 
    samtools view -bu | 
    samtools sort -o $outfile
    """

}

process qualimap{

    module 'bioinfo:qualimap/2.2.1'
    clusterOptions '--ntasks 16 --time 12:00:00 -A bharpur'

    input:
    val ngmout from ch_out_nextgenmap

    //path masterdir from masterdirchan

    output:
    val outdir into ch_out_qualimap

    script: 
    def sampID = (ngmout =~ /AGH[^\/]*/)[0]
    outdir = '/depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/qualimap/' + sampID 

    """
    echo outdir
    mkdir -p outdir

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
    sampIDstring = "\"" + qualdirs.collect { (it =~ /AGH[^\/]*/)[0] } + "\""
    qualdirString  = "\"" + qualdirs + "\""

    """
    #!/apps/spack/bell/apps/r/4.1.2-gcc-9.3.0-rw7vp7m/bin/Rscript

    sampIDstring = gsub("]","",gsub("[","",$sampIDstring,fixed=TRUE), fixed = T) 
    sampIDstring = unlist(strsplit(sampIDstring, split = ","))

    qualdirString = gsub("]","",gsub("[","",$qualdirString,fixed=TRUE), fixed = T) 
    qualdirString = unlist(strsplit(qualdirString, split = ","))

    write.table(data.frame(cbind(sampIDstring, qualdirString)),file="outframe.tsv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

    """
}

process qualimap_collate{
    
    module 'bioinfo:qualimap/2.2.1'

    input:
    file qualframe from ch_out_qualprep

    script:
    """
    mkdir -p /depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/qualimap_multiqc/
    qualimap multi-bamqc -d $qualframe -outdir /depot/bharpur/data/popgenomes/USDAHornets/flowtest/output/qualimap_multiqc/
    """
}