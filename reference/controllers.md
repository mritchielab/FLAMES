# Get controllers

Gets the controllers for the pipeline.

## Usage

``` r
controllers(pipeline)

# S4 method for class 'FLAMES.Pipeline'
controllers(pipeline)
```

## Arguments

- pipeline:

  A FLAMES.Pipeline object.

## Value

A named list of `crew_class_controller` objects, where each controller
corresponds to a step in the pipeline.

## Examples

``` r
pipeline <- example_pipeline(type = "MultiSampleSCPipeline")
#> Writing configuration parameters to:  /tmp/RtmpaNKtnp/file80d031ac0c07/config_file_32976.json 
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
controllers(pipeline) # get the controllers
#> $default
#> <crew_class_controller>
#>   Public:
#>     autoscale: function (loop = later::current_loop(), controllers = NULL) 
#>     backup: active binding
#>     cancel: function (names = character(0L), all = FALSE) 
#>     client: active binding
#>     collect: function (scale = TRUE, throttle = TRUE, error = NULL, controllers = NULL) 
#>     crashes: function (name, controllers = NULL) 
#>     crashes_max: active binding
#>     descale: function (controllers = NULL) 
#>     empty: function (controllers = NULL) 
#>     error: active binding
#>     garbage_collection: active binding
#>     initialize: function (client = NULL, launcher = NULL, reset_globals = NULL, 
#>     launch: function (n = 1L, controllers = NULL) 
#>     launcher: active binding
#>     loop: active binding
#>     map: function (command, iterate, data = list(), globals = list(), 
#>     nonempty: function (controllers = NULL) 
#>     pids: function (controllers = NULL) 
#>     pop: function (scale = TRUE, collect = NULL, throttle = TRUE, error = NULL, 
#>     pop_backlog: function (controllers = NULL) 
#>     private: environment
#>     profile: active binding
#>     push: function (command, data = list(), globals = list(), substitute = TRUE, 
#>     push_backlog: function (name, controller = NULL) 
#>     queue_backlog: active binding
#>     queue_resolved: active binding
#>     reset_globals: active binding
#>     reset_options: active binding
#>     reset_packages: active binding
#>     resolved: function (controllers = NULL) 
#>     saturated: function (collect = NULL, throttle = NULL, controller = NULL) 
#>     scale: function (throttle = TRUE, controllers = NULL) 
#>     self: crew_class_controller, R6
#>     size: function (controllers = NULL) 
#>     start: function (controllers = NULL) 
#>     started: function (controllers = NULL) 
#>     summary: function (controllers = NULL) 
#>     tasks: active binding
#>     terminate: function (controllers = NULL) 
#>     unresolved: function (controllers = NULL) 
#>     validate: function () 
#>     wait: function (mode = "all", seconds_interval = NULL, seconds_timeout = Inf, 
#>     walk: function (command, iterate, data = list(), globals = list(), 
#>   Private:
#>     .backup: NULL
#>     .client: crew_class_client, R6
#>     .crash_log: environment
#>     .crashes_max: 5
#>     .error: NULL
#>     .garbage_collection: FALSE
#>     .launcher: crew_class_launcher_local, crew_class_launcher, R6
#>     .loop: NULL
#>     .name_new_task: function (name) 
#>     .queue_backlog: NULL
#>     .queue_resolved: NULL
#>     .register_started: function () 
#>     .reset_globals: TRUE
#>     .reset_options: FALSE
#>     .reset_packages: FALSE
#>     .resolve: function (force) 
#>     .scan_crash: function (name, task) 
#>     .summary: NULL
#>     .tasks: environment
#> 
```
