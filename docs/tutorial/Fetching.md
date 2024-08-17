# Fetching public NGS Data

Let's say you want to fetch the data from the public server ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) using the ID [GSE52856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52856), which is a dataset of RNA-seq data from N2a mouse cell lines treated with shRNA against the gene *Ptbp1* and/or *Ptbp2*, and store the data in the directory `/path/to/output`. The command will be:

``` bash
# Docker
docker run --rm -v $(PWD):$(PWD) naotokubota/ngsfetch:latest \
ngsFetch -i GSE52856 -o $(PWD)/GSE52856

# Singularity
singularity exec docker://naotokubota/ngsfetch:latest \
ngsFetch -i GSE52856 -o ./GSE52856
```

`./GSE52856` will contain the data fetched from the public server:

``` bash
GSE52856/
├── fastq
│   ├── log
│   │   ├── md5sum.log
│   │   ├── SRR1043121_1.fastq.gz.log
│   │   ├── SRR1043121_2.fastq.gz.log
│   │   ├── SRR1043122_1.fastq.gz.log
│   │   ├── SRR1043122_2.fastq.gz.log
│   │   ├── SRR1107524_1.fastq.gz.log
│   │   ├── SRR1107524_2.fastq.gz.log
│   │   ├── SRR1107525_1.fastq.gz.log
│   │   └── SRR1107525_2.fastq.gz.log
│   ├── SRR1043121_1.fastq.gz
│   ├── SRR1043121_2.fastq.gz
│   ├── SRR1043122_1.fastq.gz
│   ├── SRR1043122_2.fastq.gz
│   ├── SRR1107524_1.fastq.gz
│   ├── SRR1107524_2.fastq.gz
│   ├── SRR1107525_1.fastq.gz
│   └── SRR1107525_2.fastq.gz
├── GSE52856_fastq_ftp.txt
└── GSE52856_metadata.json
```

`fastq` directory contains the fastq files fetched from the public server, and `GSE52856_fastq_ftp.txt` and `GSE52856_metadata.json` are the files containing the FTP links and metadata of the fetched data, respectively. `md5sum.log` files are the log files for the MD5 check of the fetched data. If the MD5 check succeeds, `md5sum.log` will contain the message below:

``` bash
./GSE52856/fastq/SRR1043121_1.fastq.gz: OK
./GSE52856/fastq/SRR1043121_2.fastq.gz: OK
./GSE52856/fastq/SRR1043122_1.fastq.gz: OK
./GSE52856/fastq/SRR1043122_2.fastq.gz: OK
./GSE52856/fastq/SRR1107524_1.fastq.gz: OK
./GSE52856/fastq/SRR1107524_2.fastq.gz: OK
./GSE52856/fastq/SRR1107525_1.fastq.gz: OK
./GSE52856/fastq/SRR1107525_2.fastq.gz: OK
```
