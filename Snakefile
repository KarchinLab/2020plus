from os.path import join

# configuration file
configfile: "config.yaml"

# output directory
output_dir=config["output_dir"]
# MAF file containing mutations
mutations=config["mutations"]
# pre-trained classifier
trained_classifier=config["trained_classifier"]

# params for simulations
num_iter=10
ids=list(map(str, range(1, num_iter+1)))

# minimum recurrent missense
min_recur=3

###################################
# Top-level rules
###################################

rule all:
    input: join(output_dir, "output/results/r_random_forest_prediction.txt")

# same rule is "all", but semantically more meaningful
rule predict:
    input: join(output_dir, "output/results/r_random_forest_prediction.txt")

# top-level rule to only train the 20/20+ random forest
rule train:
    input: join(output_dir, "2020plus.Rdata")

# use an already trained 20/20+ random forest to predict new data
rule pretrained_predict:
    input: join(output_dir, "pretrained_output/results/r_random_forest_prediction.txt")

###################################
# Code for calculating empirical null
# distribution based on simulations
###################################

# Simulate MAF files for subsequent running by oncogene/tsg test
rule simMaf:
    input:
        MUTATIONS=mutations
    params:
        min_recur=min_recur,
        data_dir=config["data_dir"]
    output: 
        join(output_dir, "simulated_summary/chasm_sim_maf{iter}.txt")
    shell:
        "mut_annotate --log-level=INFO "
        "  -b {params.data_dir}/snvboxGenes.bed -i {params.data_dir}/snvboxGenes.fa -c 1.5 "
        "  -m {input.MUTATIONS} -p 0 -n 1 --maf --seed=$(({wildcards.iter}*42)) "
        "  -r {params.min_recur} --unique -o {output}"

# calculate summarized features for the simulated mutations
rule simSummary:
    input: 
        MUTATIONS=mutations
    params:
        min_recur=min_recur,
        data_dir=config["data_dir"]
    output: 
        join(output_dir, "simulated_summary/chasm_sim_summary{iter}.txt")
    shell:
        "mut_annotate --log-level=INFO "
        "  -b {params.data_dir}/snvboxGenes.bed -i {params.data_dir}/snvboxGenes.fa "
        "  -c 1.5 -m {input.MUTATIONS} -p 0 -n 1 --summary --seed=$(({wildcards.iter}*42)) "
        "  --score-dir={params.data_dir}/scores "
        "  --unique -r {params.min_recur} -o {output}"


# run probabilistic2020 tsg statistical test on simulated MAF
rule simTsg:
    input:
        join(output_dir, "simulated_summary/chasm_sim_maf{iter}.txt")
    params:
        num_sim=config["NUMSIMULATIONS"],
        data_dir=config["data_dir"]
    threads: 10
    output:
        join(output_dir, "simulated_summary/tsg_sim{iter}.txt")
    shell:
        "probabilistic2020 --log-level=INFO tsg "
        "  -c 1.5 -n {params.num_sim} -b {params.data_dir}/snvboxGenes.bed "
        "  -m {input} -i {params.data_dir}/snvboxGenes.fa -p {threads} -d 1 "
        "  -o {output} "

# run probabilistic2020 oncogene statistical test on simulated MAF
rule simOg:
    input:
        mutations=join(output_dir, "simulated_summary/chasm_sim_maf{iter}.txt")
    params:
        min_recur=min_recur,
        num_sim=config["NUMSIMULATIONS"],
        data_dir=config["data_dir"]
    threads: 10
    output:
        join(output_dir, "simulated_summary/oncogene_sim{iter}.txt")
    shell:
        "probabilistic2020 --log-level=INFO oncogene "
        "  -c 1.5 -n {params.num_sim} -b {params.data_dir}/snvboxGenes.bed "
        "  -m {input.mutations} -i {params.data_dir}/snvboxGenes.fa -p {threads} " 
        "  --score-dir={params.data_dir}/scores -r {params.min_recur} "
        "  -o {output}"

# Combine the results from simOg, simTsg, and simSummary
rule simFeatures:
    input:
        summary=join(output_dir, "simulated_summary/chasm_sim_summary{iter}.txt"),
        og=join(output_dir, "simulated_summary/oncogene_sim{iter}.txt"),
        tsg=join(output_dir, "simulated_summary/tsg_sim{iter}.txt")
    params:
        data_dir=config["data_dir"]
    output:
        join(output_dir, "simulated_summary/simulated_features{iter}.txt")
    shell:
        "python `which 2020plus.py` features "
        "  -s {input.summary} --tsg-test {input.tsg} -og-test {input.og} "
        "  --permute-biogrid -o {output}"

# final processing of the simulation results
rule finishSim:
    input:
        expand(join(output_dir, "simulated_summary/simulated_features{iter}.txt"), iter=ids)
    output:
        join(output_dir, "simulated_summary/simulated_features.txt")
    shell:
        'cat {input} | awk -F"\t" \'{{OFS="\t"}} NR == 1 || !/^gene/\' - > ' + output_dir + '/simulated_summary/tmp_simulated_features.txt ; '
        'cat '+output_dir+'/simulated_summary/tmp_simulated_features.txt | awk -F"\t" \'{{OFS="\t"}}{{if(NR != 1) printf (NR"\t"); if(NR!=1) for(i=2; i<NF; i++) printf ($i"\t"); if(NR != 1) print $i; if(NR==1) print $0}}\' - > {output}'

###################################
# Code for calculating results on 
# actually observed mutations
###################################

# calculate summarized features for the observed mutations
rule summary:
    input: 
        mutations=mutations
    params:
        min_recur=min_recur,
        data_dir=config["data_dir"]
    output: 
        join(output_dir, "summary.txt")
    shell:
        "mut_annotate --log-level=INFO "
        "  -b {params.data_dir}/snvboxGenes.bed -i {params.data_dir}/snvboxGenes.fa "
        "  -c 1.5 -m {input.mutations} -p 0 -n 0 --summary "
        "  --score-dir={params.data_dir}/scores "
        "  --unique -r {params.min_recur} -o {output}"

# run probabilistic2020 tsg statistical test on MAF
rule tsg:
    input:
        mutations
    params:
        num_sim=config["NUMSIMULATIONS"],
        data_dir=config["data_dir"]
    threads: 10
    output:
        join(output_dir, "tsg.txt")
    shell:
        "probabilistic2020 -v --log-level=INFO tsg "
        "  -c 1.5 -n {params.num_sim} -b {params.data_dir}/snvboxGenes.bed "
        "  -m {input} -i {params.data_dir}/snvboxGenes.fa -p {threads} -d 1 "
        "  -o {output} "

# run probabilistic2020 oncogene statistical test on MAF
rule og:
    input:
        mutations=mutations
    params:
        min_recur=min_recur,
        num_sim=config["NUMSIMULATIONS"],
        data_dir=config["data_dir"]
    threads: 10
    output:
        join(output_dir, "oncogene.txt")
    shell:
        "probabilistic2020 -v --log-level=INFO oncogene "
        "  -c 1.5 -n {params.num_sim} -b {params.data_dir}/snvboxGenes.bed "
        "  -m {input.mutations} -i {params.data_dir}/snvboxGenes.fa -p {threads} " 
        "  --unique --score-dir={params.data_dir}/scores -r {params.min_recur} "
        "  -o {output}"

# Combine the results from og, tsg, and summary
rule features:
    input:
        summary=join(output_dir, "summary.txt"),
        og=join(output_dir, "oncogene.txt"),
        tsg=join(output_dir, "tsg.txt")
    params:
        data_dir=config["data_dir"]
    output:
        join(output_dir, "features.txt")
    shell:
        "python `which 2020plus.py` features "
        "  -s {input.summary} --tsg-test {input.tsg} -og-test {input.og} "
        "  -o {output}"

# perform prediction by random forest
# in this case the data is pan-cancer
# and so a cross-validation loop is performed
rule cv_predict:
    input:
        features=join(output_dir, "features.txt"),
        sim_features=join(output_dir, "simulated_summary/simulated_features.txt"),
    params:
        ntrees=200,
        data_dir=config["data_dir"],
        output_dir=config["output_dir"]
    output: 
        join(output_dir, "output/results/r_random_forest_prediction.txt")
    shell:
        """
        python `which 2020plus.py` --log-level=INFO train -d .7 -o 1.0 -n {{params.ntrees}} -r {outdir}/trained.Rdata --features={{input.features}} --random-seed 71
        python `which 2020plus.py` --log-level=INFO classify --trained-classifier {outdir}/trained.Rdata --null-distribution {outdir}/simulated_null_dist.txt --features {{input.sim_features}} --simulated
        python `which 2020plus.py` --out-dir {outdir}/output --log-level=INFO classify -n {{params.ntrees}} -d .7 -o 1.0 --features {{input.features}} --null-distribution {outdir}/simulated_null_dist.txt --random-seed 71
        """.format(outdir=output_dir)

#############################
# Rules for just training on
# pan-cancer data
#############################
rule train_pancan:
    input:
        features=join(output_dir, "features.txt")
    params:
        ntrees=200,
        data_dir=config["data_dir"],
        output_dir=config["output_dir"]
    output:
        join(output_dir, "2020plus.Rdata")
    shell:
        """
        python `which 2020plus.py` --log-level=INFO train -d .7 -o 1.0 -n {{params.ntrees}} --features={{input.features}} --random-seed 71 -r {outdir}/2020plus.Rdata 
        """.format(outdir=output_dir)

#############################
# Rules for predicting using
# a trained classifier on a separate
# mutation data set
#############################
rule predict_test:
    input:
        trained_classifier=trained_classifier,
        features=join(output_dir, "features.txt"),
        sim_features=join(output_dir, "simulated_summary/simulated_features.txt"),
    params:
        ntrees=200,
    output: 
        join(output_dir, "pretrained_output/results/r_random_forest_prediction.txt")
    shell:
        """
        python `which 2020plus.py` --log-level=INFO classify --trained-classifier {{input.trained_classifier}} --null-distribution {outdir}/simulated_null_dist.txt --features {{input.sim_features}} --simulated
        python `which 2020plus.py` --out-dir {outdir}/pretrained_output --log-level=INFO classify -n {{params.ntrees}} --trained-classifier {{input.trained_classifier}} -d .7 -o 1.0 --features {{input.features}} --null-distribution {outdir}/simulated_null_dist.txt --random-seed 71
        """.format(outdir=output_dir)
