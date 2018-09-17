# All content Copyright (C) 2018 Genomics plc
class WecallConfig(object):

    def __init__(self, filename):
        self.filename = filename


class ConfigFileWriter(object):

    def __init__(self, filename):
        self.__filename = filename

    def __enter__(self):
        self.__file = open(self.__filename, "w")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file.close()

    def write_config_line(self, key, value):
        self.__file.write("{} = {}\n".format(key, value))


class WecallConfigBuilder(object):

    def __init__(self, wecall_input_data, filestem):
        self.filestem = filestem
        self.__configuration = {
            "refFile": wecall_input_data.reference_filename,
            "inputs": ",".join(wecall_input_data.bam_filenames)
        }

    def with_configuration(self, key, value):
        self.__configuration[key] = value
        return self

    def build(self):
        filename = self.filestem + ".cfg"
        with ConfigFileWriter(filename) as config_writer:
            for key, value in list(self.__configuration.items()):
                config_writer.write_config_line(key=key, value=value)
        return WecallConfig(filename)
