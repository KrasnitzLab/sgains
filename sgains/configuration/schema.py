import os


def validate_existing_path(field: str, value: str, error):
    if not os.path.isabs(value):
        error(field, f"path <{value}> is not an absolute path!")
    if not os.path.exists(value):
        error(field, f"path <{value}> does not exist!")


def validate_path(field: str, value: str, error):
    if not os.path.isabs(value):
        error(field, f"path <{value}> is not an absolute path!")


sge_schema = {
    "sge_queues": {"type": "list"},
    "sge_processes": {
        "type": "integer", "default": 1,
    },
    "sge_cores": {
        "type": "integer", "default": 2,
    },
    "sge_memory": {
        "type": "string", "default": "16GB",
    },
    "sge_resource_spec": {
        "type": "string", "default": "m_mem_free=16G,h_vmem=16G",
    },
    "sge_job_extra": {"type": "list"},
}


aligner_schema = {
    "aligner_name": {"type": "string", "allowed": ["bowtie", "hisat2"]},
}


genome_schema = {
    "genome_version": {"type": "string", "allowed": ["hg19", "hg38"]},
    "genome_pristine_dir": {
        "type": "string", 
        "coerce": "abspath",
        "check_with": validate_existing_path
    },
    "genome_dir": {
        "type": "string", 
        "coerce": "abspath",
        "check_with": validate_path
    },
    "genomeindex_prefix": {
        "type": "string",
        "default": "genomeindex",
    },
}


mappable_regions_schema = {
    "mappable_read_length": {"type": "integer", "default": 100},
    "mappable_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
    },
    "mappable_file": {
        "type": "string",
        "default": "mappable_regions.txt",
    },
    "mappable_aligner_options": {
        "type": "string",
        "default": "",
    },
}


bins_schema = {
    "bins_count": {
        "type": "integer", "default": 10000,
    },
    "bins_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
    },
    "bins_file": {
        "type": "string",
        "default": "bins_boundaries.txt",
    },
}

reads_schema = {
    "reads_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
    },
    "reads_suffix": {
        "type": "string",
        "default": ".fastq.gz",
    },
}

mapping_schema = {
    "mapping_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
    },
    "mapping_suffix": {
        "type": "string",
        "default": ".rmdup.bam",
    },
    "mapping_aligner_options": {
        "type": "string",
        "default": "",
    },
}


varbin_schema = {
    "varbin_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
    },
    "varbin_suffix": {
        "type": "string",
        "default": ".varbin.txt",
    },
}


scclust_schema = {
    "scclust_case": {"type": "string"},

    "scclust_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
    },
    "scclust_cytoband_file": {
        "coerce": "abspath",
        "check_with": validate_existing_path,
    },
    "scclust_nsim": {
        "type": "integer", "default": 150,
    },
    "scclust_sharemin": {
        "type": "float", "default": 0.85,
    },
    "scclust_fdrthres": {
        "type": "integer", "default": -3,
    },
    "scclust_nshare": {
        "type": "integer", "default": 4,
    },
    "scclust_climbtoshare": {
        "type": "integer", "default": 5,
    },
}


data_10x_schema = {
    "data_10x_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
    },
    "data_10x_prefix": {
        "type": "string",
        "default": ""
    },
    "data_10x_summary": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
    },
    "data_10x_bam": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
    },
    "data_10x_bai": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
    },
}


varbin_10x_schema = {
    "varbin_10x_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
    },
    "varbin_10x_suffix": {
        "type": "string",
        "default": ".varbin.txt",
    },
}


mapping_10x_schema = {
    "mapping_10x_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
    },
    "mapping_10x_suffix": {
        "type": "string",
        "default": ".rmdup.bam",
    },
    "mapping_10x_aligner_options": {
        "type": "string",
        "default": "",
    },
}


sgains_schema = {
    "sge": {
        "type": "dict", "schema": sge_schema,
    },
    "aligner": {
        "type": "dict", "schema": aligner_schema,
    },
    "genome": {
        "type": "dict", "schema": genome_schema,
    },
    "mappable_regions": {
        "type": "dict", "schema": mappable_regions_schema,
    },
    "bins": {
        "type": "dict", "schema": bins_schema,
    },
    "reads": {
        "type": "dict", "schema": reads_schema,
    }, 
    "mapping": {
        "type": "dict", "schema": mapping_schema,
    },
    "varbin": {
        "type": "dict", "schema": varbin_schema,
    },
    "scclust": {
        "type": "dict", "schema": scclust_schema,
    },
    "data_10x": {
        "type": "dict", "schema": data_10x_schema,
    },
    "mapping_10x": {
        "type": "dict", "schema": mapping_10x_schema,
    },
    "varbin_10x": {
        "type": "dict", "schema": varbin_10x_schema,
    },
}
