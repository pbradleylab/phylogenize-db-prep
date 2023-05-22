def get_pangenome_attributes(sample, attr, pep):
    # get subsample rows that match the subsample name
    return pep.get_sample(sample)[attr]