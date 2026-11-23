def get_subsample_attributes(subsample, attribute, pep):
    subsample_rows = pep.subsample_table[pep.subsample_table["subsample"] == subsample]
    return pep.get_sample(subsample_rows["sample_name"].tolist()[0])[attribute]


def get_config_bool(value, name):
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in ["true", "false"]:
        return value.strip().lower() == "true"
    raise ValueError(name + " must be true or false")


def get_annotate_only(config, default=False):
    for section in ["mode", "annotation"]:
        if "annotate_only" in config.get(section, {}):
            return get_config_bool(config[section]["annotate_only"], section + ".annotate_only")
    return default
