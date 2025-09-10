nextflow.enable.dsl = 2

process GENERATE_TRAIN_CONFIG{
    tag "$id"
    label "process_single"
    input:
        tuple val(id), path(sdf_file), val(level)
    output:
        tuple val(id), path(sdf_file), path("config.txt"), val(level)
    script:
    def validate = params.validate ? 1 : 0
    def save_msar = params.save_msar ? 1 : 0
    """
    echo "SetClassLimit ${params.class_limit}" >> config.txt
    echo "SetIAPLimit ${params.iap_limit}" >> config.txt
    echo "BaseCreate ${level} ${id}" >> config.txt
    echo "BaseAddNewData ${sdf_file} ${params.activity}" >> config.txt
    echo "BaseTraining" >> config.txt
    if [ ${validate} ]; then
    echo "BaseValidation" >> config.txt
    fi
    if [ ${save_msar} ]; then
    echo "BaseSaveAsSAS" >> config.txt
    fi
    echo "BaseClose" >> config.txt
    """
    stub:
    """
    touch config.txt
    """
}

process TRAIN{
    tag "$id"
    container "docker.io/scottyhardy/docker-wine"
    containerOptions "--network=none -e WINEARCH=win32 -e RUN_AS_ROOT=yes"
    shell "/bin/bash"
    stageInMode 'copy'
        time {
            if (params.validate){
                if (level < 5){
                    2.hour
                }
                else if (level >= 5 && level < 7){
                    2.days
                }
                else{
                    7.days
                }
            }
            else{
                if (level < 5){
                    1.hour
                }
                else if (level >= 5 && level < 7){
                    8.hour
                }
                else{
                    2.days
                }
            }
        }
        queue {
            if (task.time <= 2.hour){
                'short'
            }
            else if (task.time > 2.hour && task.time <= 16.hour){
                'medium'
            }
            else if (task.time > 16.hour && task.time <= 6.days){
                'long'
            }
            else{
                'infinite'
            }
        }     
    input:
        tuple val(id), path(sdf_file), path(config), val(level)
        path train_program
    output:
        tuple val(id), path("*.HST"), emit: result
        tuple val(id), path(sdf_file), path("*.MSAR"), emit: model, optional: true  
    script:
    """
    WINEPREFIX=$PWD
    program_path=\$(readlink -f ${train_program} |  tr '/' '\\')
    config_path=\$(readlink -f ${config} |  tr '/' '\\')
    mkdir -p /tmp/.wine-0
    chmod -R 700 /tmp/.wine-0
    wine cmd /c "start /b /wait Z:\${program_path} Z:\${config_path}"
    exit 0
    """
    stub:
    if (params.save_msar){
        """
        touch ${id}.MSAR
        touch ${id}.HST
        """    
    }
    else{
        """
        touch ${id}.HST
        """   
    }
}

workflow{
    levels = channel.of(params.min_level..params.max_level)
    input_files = Channel.empty()
    if (params.mode == "only_train"){
        input_files = Channel.fromPath(params.input + "/*.sdf") |
                combine(levels) |
                map { file, level ->
                        def parts = file.baseName.split('_')
                        [parts[0] + "_" + parts[1] + "_" + parts[2] + "_" + parts[3] + "_" + level.toString(), file, level]
        }
        configs = GENERATE_TRAIN_CONFIG(input_files)
        train_results = TRAIN(configs, params.train_program)
    }
}