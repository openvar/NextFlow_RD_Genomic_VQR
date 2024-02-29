# Running your nextflow pipeline applet on the command line

## Run the applet help command to see the input parameters

### Run using applet name
```
dx run Nextflow_RD_GENOMIC_SC -h
```

### Run using applet ID 
As each version of an applet has a unique ID, running using the applet ID is a safer way to ensure that you are running the applet version that you intend to run

or by its applet ID:

Get applet ID:
```
 dx describe Nextflow_RD_GENOMIC_SC --json | jq -r .id
applet-GfZZ522KFxQb5KQBx0k4qQpG
```

**Help command for your applet**
```
dx run applet-GfZZ522KFxQb5KQBx0k4qQpG -h

usage: dx run applet-GfZZ522KFxQb5KQBx0k4qQpG [-iINPUT_NAME=VALUE ...]

Applet: Nextflow_RD_GENOMIC_SC

Nextflow_RD_GENOMIC_SC

Inputs:
  outdir: [-ioutdir=(string)]
        (Nextflow pipeline required) Default value:./results

  genome_file: [-igenome_file=(file)]
        (Nextflow pipeline required) Default value:genome.fasta

  genome_index_files: [-igenome_index_files=(string)]
        (Nextflow pipeline required)

  samplesheet: [-isamplesheet=(file)]
        (Nextflow pipeline required) Default value:samplesheet

         Docker Credentials: [-idocker_creds=(file)]
        Docker credentials used to obtain private docker images.

 Nextflow options
  Nextflow Run Options: [-inextflow_run_opts=(string)]
        Additional run arguments for Nextflow (e.g. -profile docker).

  Nextflow Top-level Options: [-inextflow_top_level_opts=(string)]
        Additional top-level options for Nextflow (e.g. -quiet).

  Soft Configuration File: [-inextflow_soft_confs=(file) [-inextflow_soft_confs=... [...]]]
        (Optional) One or more nextflow configuration files to be appended to the Nextflow
        pipeline configuration set

  Script Parameters File: [-inextflow_params_file=(file)]
        (Optional) A file, in YAML or JSON format, for specifying input parameter values

        Additional pipeline parameters
  Nextflow Pipeline Parameters: [-inextflow_pipeline_params=(string)]
        Additional pipeline parameters for Nextflow. Must be preceded with double dash
        characters (e.g. --foo, which can be accessed in the pipeline script using the
        params.foo identifier).

 Advanced Executable Development Options
  Debug Mode: [-idebug=(boolean, default=false)]
        Shows additional information in the job log. If true, the execution log messages
        from Nextflow will also be included.

  Resume: [-iresume=(string)]
        Unique ID of the previous session to be resumed. If 'true' or 'last' is provided
        instead of the sessionID, will resume the latest resumable session run by an applet
        with the same name in the current project in the last 6 months.

          Preserve Cache: [-ipreserve_cache=(boolean, default=false)]
        Enable storing pipeline cache and local working files to the current project. If
        true, local working files and cache files will be uploaded to the platform, so the
        current session could be resumed in the future

Outputs:
  Published files of Nextflow pipeline: [published_files (array:file)]
        Output files published by current Nextflow pipeline and uploaded to the job output
        destination.
```

## Building the run command

There are slightly different ways of specifying the input parameters depending on if they are of the `'file'` or `'string'` class *when running with the CLI*. For the UI, for parameters with class `'file'`, there is a box to click to load the file, and for parameters of class string, you will just see an empty text box where you will need to use the `'dx://project-id:/path/to/folder'` format for folder paths. 

In the output of the help command for the applet, we can see that samplesheet is of the *file* class as it has `(file)` specified in `samplesheet: [-isamplesheet=(file)]`

```
  samplesheet: [-isamplesheet=(file)]
        (Nextflow pipeline required) Default value:samplesheet
```

whereas genome_index_files is of the *string* class as it has `(string)` specified in `[-igenome_index_files=(string)]`

```
  genome_index_files: [-igenome_index_files=(string)]
        (Nextflow pipeline required)
```


If `(file)` given for parameter (i.e., the input is of the 'file' class),
 use `project-GfKZFYXKFxQz8GgXjx3xBP1P:/path/to/file` e.g
`-isamplesheet=project-GfKZFYXKFxQz8GgXjx3xBP1P:/data/sample_sheets/samplesheet_dnax_mod.tsv`

If `(string)` given for parameter (usually used for folder paths; the input is of the 'string' class),
 use `dx://project-GfKZFYXKFxQz8GgXjx3xBP1P:/path/to/folder` e.g. 
`-iraw=-igenome_index_files=dx://project-GfKZFYXKFxQz8GgXjx3xBP1P:/data/genome/*.{fasta,fasta.*}` 

For more information see the [documentation on specifying inputs and outputs to nextflow DNAnexus applets ](https://documentation.dnanexus.com/user/running-apps-and-workflows/running-nextflow-pipelines#nextflow-pipeline-executable-inputs-and-outputs)

To summarise:
 - For parameters with class string such as folder paths, *add dx:// before the project-id*
 - For parameters with class 'file' e.g. files, just use the project-id without the dx:// before it
 - For other parameters of class string  e.g a parameter named 'flowcell' (which doesn't exist in this applet) just enter the string as normal e.g. `-iflowcell='flowcell-id'`. The dx://project-id format for strings is only needed when pointing to folder paths. 


**An example run command for your applet**

[Information on running nextflow applets on dnanexus](https://documentation.dnanexus.com/user/running-apps-and-workflows/running-nextflow-pipelines#running-a-nextflow-pipeline-executable-app-or-applet)

```
dx run applet-GfZZ522KFxQb5KQBx0k4qQpG -ioutdir=./testrun -igenome_file=project-GfKZFYXKFxQz8GgXjx3xBP1P:/data/genome/Homo_sapiens_assembly38.fasta -igenome_index_files=dx://project-GfKZFYXKFxQz8GgXjx3xBP1P:/data/genome -isamplesheet=project-GfKZFYXKFxQz8GgXjx3xBP1P:/data/sample_sheets/samplesheet_dnax_mod.tsv -inextflow_run_opts='-profile docker' -ipreserve_cache=true
--destination project-GfKZFYXKFxQz8GgXjx3xBP1P:/results 
```

This will place your results in project-GfKZFYXKFxQz8GgXjx3xBP1P:/results/testrun 

## Nextflow cache on DNAnexus
As I turned on preserve_cache here, you will be able to view the work directory of the run in project-GfKZFYXKFxQz8GgXjx3xBP1P:/.nextflow_cache_db/<session-id-of-run> once the run has completed.  This is not required so you can omit this parameter if you do not want to do this as the default for this parameter is set to false if it is not specified. 

You can have up to 20 caches per project. Once the limit of 20 is reached, any new nextflow runs with preserve_cache turned on will automatically fail with an alert to delete some folders from the cache before you can run. 

In order to use the nextflow resume function on DNAnexus, you must have turned on preserve_cache in the run that you want to resume before you commenced the run. For example, if you are doing a variant calling run and anticipate that you may need to resume it in future, set -ipreserve_cache=true in the run.