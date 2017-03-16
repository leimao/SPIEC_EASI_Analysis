
# coding: utf-8

# # SPIEC-EASI Data Post-Processing

# Author: Lei Mao <br>
# Date: 3/13/2017

# Using the weights data and taxonomy annotation from [SPIEC-EASI](https://github.com/zdk123/SpiecEasi), this script is able to prepare the data for microbial network inference plot. The data for the plot contains the table of positive edge fractions and the table of edge pair numbers. The default analysis is at "Order" level. Analysis at deeper levels, such as "Family" and "Genus" may have problems due to "NA" in some OTUs.

# In[1]:

import pandas as pd


# In[2]:

weights_file = 'weights.csv'
taxonomy_file = 'tax_table.csv'
analysis_level = 'Order'


# In[3]:

def index_transform(weights):
    '''
    Transform 1-based indices in the weights dataframe to 0-based indices.
    weights: weights dataframe
    '''
    weights_transformed = weights.copy()
    weights_transformed['i'] = weights_transformed['i'].apply(lambda x: x - 1)
    weights_transformed['j'] = weights_transformed['j'].apply(lambda x: x - 1)
    return weights_transformed


# In[4]:

def taxonomy_annotate(weights, taxonomy, level):
    '''
    Add taxonomy annotations of certain level to the weights dataframe using the indices.
    weights: weights dataframe
    taxonomy: taxonomy dataframe
    level: level name
    '''
    weights_added = weights.copy()
    weights_added['tax_1'] = weights_added['i'].apply(lambda x: taxonomy.iloc[x][level])
    weights_added['tax_2'] = weights_added['j'].apply(lambda x: taxonomy.iloc[x][level])
    return weights_added


# In[5]:

def edge_count(weights, tax_1, tax_2):
    '''
    Count and calculate the positive propotions of certain taxonomy pair (tax_1 and tax_2).
    weights: weights dataframe
    tax_1: taxonomy name 1
    tax_2: taxonomy name 2
    Return: association, number of edges
    '''
    weights_interest = weights[((weights['tax_1'] == tax_1) & (weights['tax_2'] == tax_2)) | ((weights['tax_1'] == tax_2) & (weights['tax_2'] == tax_1))]
    num_positive = len(weights_interest[weights_interest['x'] > 0])
    num_negative = len(weights_interest[weights_interest['x'] < 0])
    if len(weights_interest) > 0:
        return (float(num_positive) / len(weights_interest), len(weights_interest))
    else:
        return ('NaN', 0)


# In[6]:

def network_inference(weights, taxonomy, level):
    '''
    Prepare network inference table for analysis and plot.
    weights: weights dataframe
    taxonomy: taxonomy dataframe
    level: level name
    '''
    weights_annotated = taxonomy_annotate(weights = weights, taxonomy = taxonomy, level = level)
    tax_list = taxonomy[level].unique().tolist()
    network = pd.DataFrame(index = tax_list, columns = tax_list)
    num_edge = pd.DataFrame(index = tax_list, columns = tax_list)
    for tax_1 in tax_list:
        for tax_2 in tax_list:
            count = edge_count(weights = weights_annotated, tax_1 = tax_1, tax_2= tax_2)
            network.ix[tax_1][tax_2] = count[0]
            num_edge.ix[tax_1][tax_2] = count[1]
    return network, num_edge


# In[7]:

weights = pd.read_csv(weights_file, sep=',')
weights = index_transform(weights)


# In[8]:

weights.head()


# In[9]:

tax = pd.read_csv(taxonomy_file, sep=',')


# In[10]:

tax.head()


# In[11]:

network, num_edge = network_inference(weights = weights, taxonomy = tax, level = analysis_level)


# In[12]:

network


# In[13]:

num_edge


# In[14]:

network.to_csv('network_' + analysis_level + '.csv')


# In[15]:

num_edge.to_csv('num_edge_' + analysis_level + '.csv')

