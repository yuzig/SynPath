import yaml


class Config(object):
    """Simple dict wrapper that adds a thin API allowing for slash-based retrieval of
    nested elements, e.g. cfg.get_config("meta/dataset_name")
    """
    sub_dictionary = {}

    def __init__(self, config_path):
        with open(config_path) as cf_file:
            self._data = yaml.safe_load(cf_file.read())
            self.sub_dictionary = dict(self._data)

    def get(self, path=None):
        # sub_dictionary = dict(self._data)
        return self.sub_dictionary[path]
