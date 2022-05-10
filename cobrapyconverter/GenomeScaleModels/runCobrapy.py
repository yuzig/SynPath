from GenomeScaleModels.CobraConverter import CobraConverter
import sys
import os


workingdirectory = os.getcwd()


def runner():
    model_path = sys.argv[1]
    list_of_paths = sys.argv[2]
    converter = CobraConverter(model_path)
    converter.run(list_of_paths)


runner()
