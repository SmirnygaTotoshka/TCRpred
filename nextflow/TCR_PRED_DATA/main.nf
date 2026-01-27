process FETCH{
    
    tag "${database}"
    label "process_single"
    conda "$workflow.projectDir/envs/tcr_pred.yml"

    secret 'MYSQL_DB_USER'
    secret 'MYSQL_DB_SERVER'
    secret 'MYSQL_DB_PASSWORD'

    input:
        tuple val(database), val(option)
    output:
        tuple val(database), path("*.csv")

    script:
        def optional = option == "mysql" ? "-s ${secrets.MYSQL_DB_SERVER} -u ${secrets.MYSQL_DB_USER} -p ${secrets.MYSQL_DB_PASSWORD} -d ${database}" : ""
        """
        python ${workflow.projectDir}/scripts/tcr_data_acquire.py  ${option} ${optional} -o \${PWD}
        """
    stub:
        """
        touch ${database}.csv
        """
}

process CLEAN{
    tag "${database}"
    label "process_single"
    conda "$workflow.projectDir/envs/tcr_pred.yml"

    input:
        tuple val(database), path(raw_table)
    output:
        tuple val("tcr_epitope"), path("*tcr_epi_cleaned.csv"), emit: epitope
        tuple val("tcr_mhc"), path("*tcr_mhc_cleaned.csv"), emit: mhc

    script:
        """
        python ${workflow.projectDir}/scripts/tcr_data_clean.py -d ${database} -i ${raw_table[0]} -o \${PWD}
        """
    stub:
        """
        touch ${database}_tcr_epi_cleaned.csv
        touch ${database}_tcr_mhc_cleaned.csv
        """
}

process MERGE{

    tag "${prefix}"
    label "process_single"
    conda "$workflow.projectDir/envs/tcr_pred.yml"

    input:
        tuple val(prefix), path(list_databases)
    output:
        tuple val(prefix), path("${prefix}.csv")
    script:
        """
        #!/usr/bin/env python
        import pandas as pd
        paths = "${list_databases}".split(" ")
        raw_merged = pd.concat([pd.read_csv(p, sep = ";", header = 0) for p in paths], ignore_index = True)
        merged_dedup = raw_merged.drop_duplicates(subset = ["Chain","Species","Structure","Activity"], ignore_index = True, keep = "first")
        merged_dedup.to_csv("${prefix}.csv",sep=";",header = True, index = False)
        """
    stub:
        """
        touch ${prefix}.csv
        """

}

process CALC_STAT{
    tag "${table}"
    conda "$workflow.projectDir/envs/tcr_pred.yml"

    input:
        tuple val(prefix), path(table)
    output:
        path("*.db")
    script:
    def table_name = table.baseName
    def database = table_name == prefix ? table_name : table_name.split("_")[0]
    if (prefix == "tcr_epitope"){
        """
        python ${workflow.projectDir}/scripts/tcr_data_statistics.py \
        --input ${table[0]} \
        --output ${database}_epi.db \
        --epitope \
        --threads ${task.cpus}
        
        databases=\$(awk -F ";" '{if (NR>1) print \$2}' ${table[0]} | uniq)
        for d in \$databases
        do
        head -n 1 ${table[0]} > tmp.csv
        awk -F ";" -v d="\$d" '{if (\$2 == d) {print \$0}}' ${table[0]} >> tmp.csv
        
        python ${workflow.projectDir}/scripts/tcr_data_statistics.py \
        --input tmp.csv \
        --output \${d}_epi.db \
        --epitope \
        --threads ${task.cpus}
        done
        """
    }
    else if (prefix == "tcr_mhc"){
        """
        python ${workflow.projectDir}/scripts/tcr_data_statistics.py \
        --input ${table[0]} \
        --output ${database}_mhc.db \
        --mhc \
        --threads ${task.cpus}
        
        databases=\$(awk -F ";" '{if (NR>1) print \$2}' ${table[0]} | uniq)
        for d in \$databases
        do
        head -n 1 ${table[0]} > tmp.csv
        awk -F ";" -v d="\$d" '{if (\$2 == d) {print \$0}}' ${table[0]} >> tmp.csv
        
        python ${workflow.projectDir}/scripts/tcr_data_statistics.py \
        --input tmp.csv \
        --output \${d}_mhc.db \
        --mhc \
        --threads ${task.cpus}
        done
        """
    }
    stub:
    def table_name = table.baseName
    def database = table_name == prefix ? table_name : table_name.split("_")[0]

    """
    touch ${database}_epi.db
    touch ${database}_mhc.db
    """
}

process SPLIT_SAMPLE{
    
    label "process_single"
    conda "$workflow.projectDir/envs/tcr_pred.yml"

    input:
        tuple val(prefix), path(total_table)
        each parameters
    output:
        path "*.csv"
    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    
    table = pd.read_csv("${total_table[0]}", sep = ";", header = 0)
    target_species = "${parameters[0]}"
    target_charged = "${parameters[1]}"    
    target_chain = "${parameters[2]}"
    target_activity = "${prefix}".split("_")[1]
    sample = table.query("Species == @target_species and Chain == @target_chain")
    sample.to_csv(f"{target_species}--cdr3{target_chain}-{target_charged}-{target_activity}.csv", sep = ";", index = False, header = True)
    """
    stub:
    """
    #!/usr/bin/env python
    target_species="${parameters[0]}"
    target_activity = "${prefix}".split("_")[1]
    target_chain="${parameters[2]}"
    target_charged="${parameters[1]}"
    filename=f"{target_species}-cdr3{target_chain}-{target_charged}-{target_activity}.csv"
    file = open(filename,"w")
    file.close()
    """
}


workflow {
    csv_records = Channel.fromPath(params.input) | splitCsv(sep: ";", header: true)
    datasets = csv_records | FETCH | CLEAN
        
    list_datasets = datasets.epitope.mix(datasets.mhc) | groupTuple
    merged_datasets = MERGE(list_datasets)
        
    statistics = CALC_STAT(merged_datasets)

    species = Channel.of("human", "mouse")
    charged = Channel.of("ch", "u") 
    chains = Channel.of("alpha", "beta")

    split_parameters = species.combine(charged).combine(chains).unique()
    splitted_datasets = SPLIT_SAMPLE(merged_datasets, split_parameters)
        
    publish:
        statistics = statistics
        total = merged_datasets
        splitted = splitted_datasets
}

output{
    statistics {
        path 'statistics'
    }
    
    total {
        path 'total'
    }
    
    splitted {
        path 'splitted'
    }
}