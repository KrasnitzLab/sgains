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
    "aligner_name": {
        "type": "string",
        "allowed": ["bowtie", "hisat2", "bwa"],
        "default": "bowtie",
        "meta": {
            "help": "aligner to use in sGAINS subcommands"
        },

    },
}


genome_schema = {
    "genome_version": {
        "type": "string",
        "allowed": ["hg19", "hg38"],
        "default": "hg19",
        "meta": {
            "help": "version of reference genome to use"
        },
    },
    "genome_pristine_dir": {
        "type": "string", 
        "coerce": "abspath",
        "check_with": validate_existing_path,
        "meta": {
            "help": "directory where clean copy of reference genome "
            "is located"
        },
    },
    "genome_dir": {
        "type": "string", 
        "coerce": "abspath",
        "check_with": validate_path,
        "meta": {
            "help": "genome index working directory"
        },
    },
    "genomeindex_prefix": {
        "type": "string",
        "default": "genomeindex",
        "meta": {
            "help": "genome index prefix"
        },
    },
}


mappable_regions_schema = {
    "mappable_read_length": {
        "type": "integer",
        "default": 100,
        "meta": {
            "help": "read length to use for generation of mappable regions"
        },
    },
    "mappable_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
        "meta": {
            "help": "directory where mappable regions working files are "
            "stored"
        },
    },
    "mappable_file": {
        "type": "string",
        "default": "mappable_regions.txt",
        "meta": {
            "help": "filename for mappable regions results"
        },
    },
    "mappable_aligner_options": {
        "type": "string",
        "default": "",
        "meta": {
            "help": "additional aligner options for use when computing "
            "uniquely mappable regions"
        },
    },
}


bins_schema = {
    "bins_count": {
        "type": "integer",
        "default": 10000,
        "meta": {
            "help": "number of bins"
        },
    },
    "bins_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
        "meta": {
            "help": "bins working directory"
        },
    },
    "bins_file": {
        "type": "string",
        "default": "bins_boundaries.txt",
        "meta": {
            "help": "bins boundaries filename"
        },
    },
}

reads_schema = {
    "reads_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
        "meta": {
            "help": "data directory where sequencing reads are located"
        },
    },
    "reads_suffix": {
        "type": "string",
        "default": ".fastq.gz",
        "meta": {
            "help": "reads files suffix pattern"
        },
    },
}

mapping_schema = {
    "mapping_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
        "meta": {
            "help": "data directory where mapping files are located"
        },
    },
    "mapping_suffix": {
        "type": "string",
        "default": ".rmdup.bam",
        "meta": {
            "help": "mapping files suffix pattern"
        },
    },
    "mapping_aligner_options": {
        "type": "string",
        "default": "",
        "meta": {
            "help": "additional aligner mapping options"
        },
    },
}


varbin_schema = {
    "varbin_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
        "meta": {
            "help": "varbin working directory"
        },
    },
    "varbin_suffix": {
        "type": "string",
        "default": ".varbin.txt",
        "meta": {
            "help": "varbin files suffix pattern"
        },
    },
}


scclust_schema = {
    "scclust_case": {
        "type": "string",
        "meta": {
            "help": "SCclust case name"
        },
    },

    "scclust_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_path,
        "meta": {
            "help": "SCclust working directory"
        },
    },
    "scclust_cytoband_file": {
        "coerce": "abspath",
        "check_with": validate_existing_path,
        "meta": {
            "help": "location of cyto band description file"
        },
    },
    "scclust_nsim": {
        "type": "integer", "default": 150,
        "meta": {
            "help": "SCclust number of simulations"
        },
    },
    "scclust_sharemin": {
        "type": "float", "default": 0.85,
        "meta": {
            "help": "SCclust sharemin parameter"
        },
    },
    "scclust_fdrthres": {
        "type": "integer", "default": -3,
        "meta": {
            "help": "SCclust fdrthres parameter"
        },
    },
    "scclust_nshare": {
        "type": "integer", "default": 4,
        "meta": {
            "help": "SCclust nshare parameter"
        },
    },
    "scclust_climbtoshare": {
        "type": "integer", "default": 5,
        "meta": {
            "help": "SCclust climbtoshare parameter"
        },
    },
}


data_10x_schema = {
    "data_10x_dir": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
        "meta": {
            "help": "10x Genomics dataset directory"
        },
    },
    "data_10x_prefix": {
        "type": "string",
        "default": "",
        "meta": {
            "help": "10x Genomics common prefix"
        },
    },
    "data_10x_cell_summary": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
        "meta": {
            "help": "10x Genomics per cell summary filename"
        },
    },
    "data_10x_bam": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
        "meta": {
            "help": "10x Genomics dataset aligned reads data"
        },
    },
    "data_10x_bai": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
        "meta": {
            "help": "10x Genomics dataset aligned reads index file"
        },
    },
}


# varbin_10x_schema = {
#     "varbin_10x_dir": {
#         "type": "string",
#         "coerce": "abspath",
#         "check_with": validate_path,
#     },
#     "varbin_10x_suffix": {
#         "type": "string",
#         "default": ".varbin.txt",
#     },
# }

# reads_10x_schema = {
#     "reads_10x_dir": {
#         "type": "string",
#         "coerce": "abspath",
#         "check_with": validate_path,
#         "meta": {
#             "help": "data directory where sequencing reads are located"
#         },
#     },
#     "reads_10x_suffix": {
#         "type": "string",
#         "default": ".fastq.gz",
#         "meta": {
#             "help": "reads files suffix pattern"
#         },
#     },

# }

# mapping_10x_schema = {
#     "mapping_10x_dir": {
#         "type": "string",
#         "coerce": "abspath",
#         "check_with": validate_path,
#     },
#     "mapping_10x_suffix": {
#         "type": "string",
#         "default": ".rmdup.bam",
#     },
#     "mapping_10x_aligner_options": {
#         "type": "string",
#         "default": "",
#     },
# }


sgains_schema = {
    "work_dirname": {
        "type": "string",
        "coerce": "abspath",
        "check_with": validate_existing_path,
    },
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
    # "reads_10x": {
    #     "type": "dict", "schema": reads_10x_schema,
    # },
    # "mapping_10x": {
    #     "type": "dict", "schema": mapping_10x_schema,
    # },
    # "varbin_10x": {
    #     "type": "dict", "schema": varbin_10x_schema,
    # },
}
