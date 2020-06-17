# Get MAFs for Neanderthal-introgressed SNPS from 1000 genomes all populations data
This document describes the process for extracting minor allele frequencies (MAFs) for Neanderthal-introgressed SNPs from the 1000 genomes all populations dataset, which comprises [26 populations](https://www.internationalgenome.org/category/population/).

Neanderthal SNPs from:
1. Dannemann M, Prufer K & Kelso J. Functional implications of Neandertal introgression in modern humans. Genome Biol 2017 18:61.
2. Simonti CN et al. The phenotypic legacy of admixture between modern humans and Neandertals. Science 2016 351:737-41.

1000 genomes data from:
* [1000 genomes](https://www.internationalgenome.org/data/)

The combined list of Neanderthal-introgressed SNPs from both sources is saved in .csv format under: /well/jknight/shiyao/data/comparison_df.csv

First, a script was written to define a "get_pop" function that extracts MAFs for Neanderthal-introgressed SNPs present in the 1000 genomes dataset from one population.

```python
def get_pop(pop='YRI'):
    import pandas as pd
    
    # Load Neanderthal-introgressed SNPs dataset
    comparison_df = pd.read_csv("/well/jknight/shiyao/data/comparison_df.csv")
    comparison_df = comparison_df.drop('Unnamed: 0', axis=1)
    comparison_df['Source'] = comparison_df['Source'].astype('str')
    pop_df = comparison_df

    # Loop over chromosomes 1-22
    starter = True
    for i in range(1, 23):
        file_name = '/well/jknight/shiyao/data/1000genome/all_populations/' + pop + '/' + pop + '_chr' + str(i) + '_freq_genome.txt'
        
        # Load population dataset
        df = pd.read_csv(file_name, sep='\t', header=None)
        df.drop([4, 5], axis=1, inplace=True)
        df = df.rename(columns={0: "Location", 1: "ID", 2: "Major", 3: "Minor", 6: pop})
        
        # Split Location column into Chromosome & Position columns
        new = df['Location'].str.split("-", expand=True)
        df['Chromosome'] = new[0]
        df["Position"] = new[1]
        df.drop("Location", axis=1, inplace=True)
        
        # Format df
        df = df[['Chromosome', 'Position', 'ID', 'Major', 'Minor', pop]]
        df[['Chromosome','Position']] = df[['Chromosome','Position']].astype('int64')
        df[['ID', 'Major', 'Minor']] = df[['ID', 'Major', 'Minor']].astype('str')  # 'string' in pandas 1.0.3
        # Sort columns with multiple values
        df.Minor = df.Minor.str.split(',').apply(sorted, 1).str.join(',').str.strip(',')
        df.ID = df.ID.str.split(';').apply(sorted, 1).str.join(';').str.strip(';')

        # Merge dataframes
        if starter:
            pop_df = pop_df.merge(df, how='left', on=['Chromosome', 'Position'])
            starter = False
        else:
            pop_df = pop_df.merge(df, how='left', on=['Chromosome', 'Position'])
            
            # Fill NA values
            pop_df['ID_x'] = pop_df['ID_x'].fillna(pop_df['ID_y'])
            pop_df['Major_x'] = pop_df['Major_x'].fillna(pop_df['Major_y'])
            pop_df['Minor_x'] = pop_df['Minor_x'].fillna(pop_df['Minor_y'])
            pop_df[pop + '_x'] = pop_df[pop + '_x'].fillna(pop_df[pop + '_y'])
            pop_df.drop(['ID_y', 'Major_y', 'Minor_y', pop + '_y'], inplace=True, axis=1)
            
            pop_df.rename(columns={'ID_x': 'ID', 'Major_x': 'Major', 'Minor_x': 'Minor', pop + '_x': pop}, inplace=True)

    # Save pop df
    pop_df.to_csv("/well/jknight/shiyao/data/1000genome/pop_df/pop_"+pop+"_df.csv")
 ```
 
Next, the following script was executed to create a *.sh file to run the "get_pop" function for each of the 26 populations and submit them to the cluster as separate jobs.

```python
import subprocess
from time import sleep
from glob import glob

# Get list of populations
populations = []
files = glob('/well/jknight/shiyao/data/1000genome/all_populations/*')
for file_name in files:
    temp = file_name.split('/')
    populations.append(temp[-1])
    
# Loop over list of populations
for pop in populations:

    # Write the *.sh file to be used in the job
    f = open(pop + '.sh', "w")

    f.write('#!/bin/bash\n\n')

    f.write('#$ -cwd -V\n')
    f.write('#$ -N ' + pop + ' -j y\n')
    f.write('#$ -P jknight.prjc -q short.qc\n')
    f.write('#$ -pe shmem 1\n')
    f.write('#$ -o ' + pop + '.out -e ' + pop + '.err\n')
    f.write('#$ -r y\n\n')

    f.write('echo "#######################################################################"\n')
    f.write('echo "SGE job id: "$JOB_ID\n')
    f.write('echo "Run on host: "`hostname`\n')
    f.write('echo "Operating system: "`uname -s`\n')
    f.write('echo "Username: "`whoami`\n')
    f.write('echo "Started at: "`date`\n')
    f.write('echo "#######################################################################"\n\n')

    # Run get_pop function
    f.write("""/apps/well/python/3.5.2-gcc5.4.0/bin/./python -c "import getpop; getpop.get_pop(pop='"""+pop+"""')"\n\n""")

    f.write('echo "#########################################################################"\n')
    f.write('echo "Finished at: "`date`\n')
    f.write('echo "#########################################################################"\n')
    f.write('exit 0')

    f.close()

    # Submit the job and grab the job_id for use in making the next job a dependency
    submit_name = 'qsub ' + pop + '.sh'
    job = subprocess.Popen([submit_name], shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,)
    output, err = job.communicate()
    print(output)
    sleep(0.5)
```

Finally, the individual pop_df.csv files were combined into a single allpop_df.csv file.

```python
from glob import glob
import pandas as pd

# Get list of populations
populations = []
files = glob('/well/jknight/shiyao/data/1000genome/all_populations/*')
for file_name in files:
    temp = file_name.split('/')
    populations.append(temp[-1])
populations.remove('CHS') # CHS population data corrupted

# Create list of pop_dfs
dfs = []
for pop in populations:
    df = pd.read_csv("/well/jknight/shiyao/data/1000genome/pop_df/pop_" + pop + "_df.csv")
    df = df.drop('Unnamed: 0', axis=1)
    df[['Source', 'ID', 'Major', 'Minor']] = df[['Source', 'ID', 'Major', 'Minor']].astype('str')
    dfs.append(df)

# Extract MAF columns from each population and add to allpop_df
starter = True
for df in dfs:
    if starter:
        allpop_df = df
        starter = False
    else:
        allpop_df = allpop_df.merge(df, how='left', on=['Chromosome', 'Position', 'Source', 'ID', 'Major', 'Minor'])

# Save allpop_df to csv
allpop_df.to_csv("/well/jknight/shiyao/data/allpop_df.csv")
```
