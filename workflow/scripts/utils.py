import math


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


def get_linclust_settings(config):
    """Validate clustering criteria, preserving the existing config key spelling."""
    options = config.get("mmseqs2", {}).get("linclust", {})
    settings = {}
    for key, name, default in [
        ("seq_id_precent", "min_seq_id", 0.0),
        ("coverage", "coverage", 0.8),
        ("evalue", "evalue", 0.001),
    ]:
        value = options.get(key, default)
        valid = isinstance(value, (int, float)) and not isinstance(value, bool)
        if valid:
            valid = math.isfinite(value) and value >= 0
            if key != "evalue":
                valid = valid and value <= 1
        if not valid:
            bounds = "a finite nonnegative number" if key == "evalue" else "a number between 0 and 1"
            raise ValueError(f"mmseqs2.linclust.{key} must be {bounds}")
        settings[name] = value

    cov_mode = options.get("cov_mode", 0)
    if type(cov_mode) is not int or cov_mode not in range(6):
        raise ValueError("mmseqs2.linclust.cov_mode must be an integer from 0 to 5")
    settings["cov_mode"] = cov_mode
    return settings
