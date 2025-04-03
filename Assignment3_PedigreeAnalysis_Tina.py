"""
Assignment 3: Pedigree Analysis
Authors: Patrick McGrath, Tucker Lancaster

In this assignment, you will practice using scikit-allel and numpy to analyze genetic variant data in order to
identify the sex-determining locus in a species of cichlid. Your task is to fill in the missing code in sections 2, 3,
and 4, and to answer the open-ended question in section 5. Missing code will be either in the form of a variable
set equal to some placeholder (like None or []) that you need to replace, or a comment that says "# YOUR CODE HERE".
"""

"""
Section 1: Setup

This section is already complete, but make sure you understand it before moving on as many of these variables will
be reused. If you are missing any of the libraries, be sure to install them.
"""

import allel, pdb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


vcf_path = './YHPedigree_Class.vcf'    # modify this line if needed to match your file location
vcf_dict = allel.read_vcf(vcf_path, fields=['variants/CHROM', 'calldata/GT', 'calldata/DP', 'samples'])
chromosome = vcf_dict['variants/CHROM']
genotype = vcf_dict['calldata/GT']
depth = vcf_dict['calldata/DP']
samples = vcf_dict['samples']


"""
Section 2: Mendelian Phasing

This section will guide you through writing code for performing haplotype phasing of the provided VCF data
"""
# First, convert the "genotype" object from a numpy.ndarray to an allel.GenotypeArray object.
# Replace "None" with your code
converted_genotype = allel.GenotypeArray(genotype)

# Next, use array indexing to isolate the paternal and maternal genotypes from the larger "converted_genotypes" array.
# Remember that these correspond to the first and second columns of the array, respectively.
# Replace "None" with your code
paternal_genotype = converted_genotype[:, 0]
maternal_genotype = converted_genotype[:, 1]

# Create a 1d boolean array that is True where the paternal genotype is heterozygous AND the maternal genotype is
# homozygous, and False otherwise. Hint: use the is_het() and is_hom() methods of the allel.GenotypeArray class
# and the & operator (& is shorthand for the np.logical_and() function)
# Replace "None" with your code
paternal_mask = paternal_genotype.is_het() & maternal_genotype.is_hom()

# create another 1d boolean array that is instead True where the maternal genotype is heterozygous AND the paternal
# genotype is homozygous, and False otherwise.
# Replace "None" with your code
maternal_mask = maternal_genotype.is_het() & paternal_genotype.is_hom()

# Use the maternal_mask and paternal_mask arrays to extract two sub-arrays from the larger converted_genotype array.
# These new arrays should have the same number of columns as converted_genotype, but only contain the rows where the
# corresponding mask array was True. If you are unfamiliar with boolean array indexing in numpy, see this page of the
# documentation: https://numpy.org/doc/stable/user/basics.indexing.html#boolean-array-indexing.
# Replace "None" with your code
paternal_masked_genotype = converted_genotype[paternal_mask]
maternal_masked_genotype = converted_genotype[maternal_mask]

# phase the masked genotypes using the allel.phase_by_transmission function with the window_size set to 100
# Replace "None" with your code
paternal_phased = allel.phase_by_transmission(paternal_masked_genotype, window_size=100)
maternal_phased = allel.phase_by_transmission(maternal_masked_genotype, window_size=100)

# Convert the phased maternal and paternal data to allel.HaplotypeArray objects
# Hint: use the .to_haplotypes() method of the allel.GenotypeArray objects
# Replace "None" with your code
paternal_haplotypes = paternal_phased.to_haplotypes()
maternal_haplotypes = maternal_phased.to_haplotypes()

# Uncomment the next three lines of completed code to paint offspring haplotypes by parental haplotypes.
num_progeny = samples.size - 2
paternally_painted = allel.paint_transmission(paternal_haplotypes[:,(0,1)],paternal_haplotypes[:,range(4,num_progeny*2+4,2)])
maternally_painted = allel.paint_transmission(maternal_haplotypes[:,(2,3)],maternal_haplotypes[:,range(5,num_progeny*2+4,2)])

"""
Section 3: Additional filtering

At this point, the paternally_painted and maternally_painted arrays should contain integers ranging from 1 to 7. 
(see the allel.paint_transmission documentation for more details). In this section, you will write code to further
filter these arrays. 
"""

#### first, subtract 1 from each value in the paternally_painted and maternally_painted arrays so that the values of
# interest become 0 and 1
# Replace "None" with your code
paternally_painted -= 1
maternally_painted -= 1

# next, convert the dtype of  both arrays to float so that we can have nan values
# Replace "None" with your code
paternally_painted = paternally_painted.astype(float)
maternally_painted = maternally_painted.astype(float)

# using boolean indexing, set all values greater than 1 to np.nan
# YOUR CODE HERE
paternally_painted[paternally_painted > 1] = np.nan
maternally_painted[maternally_painted > 1] = np.nan

# using the paternal_phased and maternal_phased arrays, and the .is_phased attribute, set all non-phased values in the
# paternally_painted and maternally_painted arrays to np.nan. Hint: you can use the ~ operator, or the np.logical_not
# function, to invert a boolean array.
# YOUR CODE HERE
paternally_painted[~paternal_phased.is_phased[:,2:]] = np.nan
maternally_painted[~maternal_phased.is_phased[:,2:]] = np.nan

# set all values with a sequencing depth less than 5 in the paternally_painted and maternally_painted arrays to np.nan
# HINT: you will need the array called "depth" that was created in section 1, as well as the paternal_mask
# and maternal_mask arrays from earlier in this section.
# YOUR CODE HERE

# paternally_painted[depth[paternal_mask] < 5] = np.nan
# maternally_painted[depth[maternal_mask] < 5] = np.nan
paternal_depth = depth[paternal_mask]
paternally_painted[paternal_depth[:,2:]< 5] = np.nan

maternal_depth = depth[maternal_mask]
maternally_painted[maternal_depth[:,2:] < 5] = np.nan

"""
Section 4: Prep for Plotting

In this section, you will calculate summary statistics from your painted arrays, which will then be plotted in section 
5.
"""

# first, find a list (or other iterable) of the unique chromosome names. There should be 22. Replace the [] below with
# appropriate code
unique_chromosomes = np.unique(chromosome)

# The next two lines of code do not need to be modified. They simply set up empty lists that you'll use later
all_paternal_means = []
all_maternal_means = []

# # Now fill in the missing sections of this loop to calculate summary statistics for each chromosome:
# for chrom in unique_chromosomes:
#     # extract the rows of the paternally_painted array associated with "chrom"
#     # Replace "None" with your code
#     paternal_slice = paternally_painted[chromosome == chrom]
#     # find the mean of this array along the 0th dimension, ignoring nan values. The resulting array should have 24
#     # values, each representing the summary statistic for a single individual
#     # Replace "None" with your code
#     pat_means = np.nanmean(paternal_slice, axis=0)
#     # append the pat_mean array to the all_paternal_means list
#     # YOUR CODE HERE
#     all_paternal_means.append(pat_means)
#
#     # repeat the above steps (extract, summarize, and append) with the maternally_painted array. Be sure to append to
#     # the all_maternal_means list, and to use the appropriate maternal equivalents when performing the first filtering
#     # step
#     # YOUR CODE HERE
#     maternal_slice = maternally_painted[chromosome == chrom]
#     mat_means = np.nanmean(maternal_slice, axis=0)
#     all_maternal_means.append(mat_means)

for chrom in unique_chromosomes:
    # Get indices for current chromosome
    chrom_idx = chromosome == chrom

    # Paternal means
    pat_chrom_mask = chrom_idx & paternal_mask
    pat_slice = paternally_painted[pat_chrom_mask[paternal_mask]]
    pat_means = np.nanmean(pat_slice, axis=0)
    all_paternal_means.append(pat_means)

    # Maternal means
    mat_chrom_mask = chrom_idx & maternal_mask
    mat_slice = maternally_painted[mat_chrom_mask[maternal_mask]]
    mat_means = np.nanmean(mat_slice, axis=0)
    all_maternal_means.append(mat_means)

"""
Section 5: Plotting

Once you have finished sections 2, 3, and 4, uncomment the below code to visualize your results. You do not need
to add any code to this section, but please answer the following question based on the visualizations produced (look for
a new file called assignment3_visualization.pdf after running this section.) Limit your answer to 1-2 sentences.

Open ended question: Based on the visualization, which linkage group (e.g., LG1, LG2, etc.) likely contains the
sex-determining locus, and from which parent is it inherited? How do you know? 

Your answer: LG10 looks like contains the sex-determining locus. Based on the visualization, it is inherited from the first 
parent, which is the maternal side. Because it shows distinct inheritance pattern based on sex and the most distinct 
variation we have found between different sexes are in LG10, in the first graph. Thus we can conclude that it is 
the one that contains sex-determining locus and it is inherited from maternal side. Based on the example on Canvas, it
is the paternal side. 

"""

linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4',
                 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 'NC_036786.1':'LG7', 'NC_036787.1':'LG8',
                 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11', 'NC_036791.1':'LG12',
                 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16',
                 'NC_036796.1':'LG17', 'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20',
                 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

male_ids = ['YH_025', 'YH_026', 'YH_027', 'YH_028', 'YH_029', 'YH_030', 'YH_037', 'YH_038', 'YH_039', 'YH_040', 'YH_041', 'YH_042']
female_ids = ['YH_016', 'YH_017', 'YH_018', 'YH_019', 'YH_020', 'YH_021', 'YH_022', 'YH_023', 'YH_024', 'YH_031', 'YH_032', 'YH_033']

def convert_means_list_to_df(all_means_list):
    df = pd.DataFrame(np.vstack(all_means_list), columns=samples[2:])
    df['linkage_group'] = [linkageGroups[x] for x in unique_chromosomes]
    df = df.melt(id_vars='linkage_group')
    df['sex'] = df.variable.apply(
        lambda x: 'male' if x in male_ids else 'female' if x in female_ids else 'unknown')
    df = df.drop(columns=['variable'])
    df.rename(columns={'value': 'match_quality'}, inplace=True)
    return df

paternal_df = convert_means_list_to_df(all_paternal_means)
maternal_df = convert_means_list_to_df(all_maternal_means)

fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
palette = sns.color_palette(['#89CFF0', '#F4C2C2'])
sns.boxplot(paternal_df, x='linkage_group', y='match_quality', hue='sex', ax=axes[0], palette=palette)
axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
axes[0].set(title='paternal inheritance')
sns.boxplot(maternal_df, x='linkage_group', y='match_quality', hue='sex', ax=axes[1], palette=palette)
axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
axes[0].set(title='maternal inheritance')
fig.tight_layout()
fig.savefig('assignment3_visualization.pdf')
plt.close(fig)
