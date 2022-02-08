#!/usr/bin/env nextflow

// Configurable variables
params.samplelist = '/depot/bharpur/data/popgenomes/USDAHornets/trimmomatic/hornet_IDs.txt'
params.rawdatadir = '/depot/bharpur/data/popgenomes/USDAHornets/raw_data'
params.genome = '/depot/bharpur/data/ref_genomes/VMAN/GCF_014083535.2_V.mandarinia_Nanaimo_p1.0_genomic.fna'

samplelist = params.samplelist   
rawdatadir = params.rawdatadir   
genome = params.genome

// Validate inputs
if ( samplelist.isEmpty () ) {
    exit 1, 'Please specify --samplelist with a list of expected sample IDS'
} else if ( rawdatadir.isEmpty () ) {
    exit 1, 'Please specify --rawdatadir with the directory containing raw data directories'
}

//pull the list of samples and assign to an array
File samplels = new File(samplelist)
def samples = samplels.readLines()

//for testing purposes, reassign samples as a single sample
samples = [samples.get(0), samples.get(1)]
println samples

//check that each sample has a corresponding directory in master
rawdatadirs = "ls $rawdatadir".execute().text
def mastercontains = { it -> rawdatadirs.contains(it) }
if ( !samples.every(mastercontains) ) {
    //println "Every sample in samplelist should have a corresponding folder in the master directory."
    exit 1, 'Every sample in samplelist should have a corresponding folder in the master directory.'
}

//check that each sample directory contains exactly two fastq files
def countfastas = { it -> new File("$rawdatadir/$it").listFiles() .findAll { it.name =~ /(?i)fq|fastq/ } .size() }
if ( !samples.every(it -> countfastas(it) == 2) ) {
    //println "Every sample directory should contain exactly two fastq files."
    exit 1, 'Every sample directory should contain exactly two fastq files.'
}

samples_in_trimmomatic = Channel.fromList(samples)
rawdatadirchan = Channel.value( rawdatadir )
genomechan = Channel.value( genome )

process trimmomatic {
    
    module 'bioinfo:trimmomatic/0.39'

    input:
    val sampleID from samples_in_trimmomatic
    val rawpath from rawdatadirchan

    output:
    tuple val(outid), val(sampleID) into ch_out_trimmomatic

    script:

    fq_1 = rawpath + '/' + sampleID + '/*_1.fq.gz'
    fq_2 = rawpath + '/' + sampleID + '/*_2.fq.gz'
    outdir = 'output/trimmomatic/' + sampleID
    outid = outdir + '/' + sampleID + ".fq.gz"

    println fq_1

    //TODO: re-parameterize the trimmomatic command for optimal trimming
    """
    mkdir -p $outdir

    trimmomatic \
    PE \
    $fq_1 \
    $fq_2 \
    -baseout $outid \
    LEADING:3 TRAILING:3 MAXINFO:36:0.7 MINLEN:36
    """
}

process nextgenmap{

    module 'bioinfo:samtools/1.12'
    memory '10 GB'
    clusterOptions '--ntasks 16 --time 30:00 -A bharpur'

    input:
    tuple val(trimpath), val(sampleID) from ch_out_trimmomatic
    path rawdatadir from rawdatadirchan
    path genome from genomechan

    output: 
    tuple val(outfile), val(sampleID) into ch_out_nextgenmap

    script:
    //def sampID = (trimpath =~ /AGH[^\/]*/)[0]

    trimfile1 = trimpath.replaceAll(".fq", "_1P.fq")
    trimfile2 = trimpath.replaceAll(".fq", "_2P.fq")

    outdir = 'output/nextgenmap/' + sampleID
    outfile = outdir + "/" + sampleID + "_ngm_sorted.bam"

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
    tuple val(ngmout), val(sampleID) from ch_out_nextgenmap

    output:
    tuple val(outdir), val(sampleID) into ch_out_qualimap

    script: 
    //def sampID = (ngmout =~ /AGH[^\/]*/)[0]
    outdir = 'output/qualimap/' + sampleID 

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
    tuple val(qualdirs), val(sampleIDs) from ch_out_qualimap.collect()

    output:
    file 'outframe.tsv' into ch_out_qualprep

    script:
    //sampIDstring = "\"" + qualdirs.collect { (it =~ /AGH[^\/]*/)[0] } + "\""
    sampIDstring = "\"" + sampleIDs + "\""
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
    mkdir -p output/qualimap_multiqc/
    qualimap multi-bamqc -d $qualframe -outdir output/qualimap_multiqc/
    """
}