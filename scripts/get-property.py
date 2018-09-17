# All content Copyright (C) 2018 Genomics plc
import json
import sys


def main(filename, property):
    with open(filename, 'r') as fp:
        props = json.load(fp)
    print(props[property])


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise Exception(
            "Usage:\n\t{exe} filename property".format(
                exe=sys.argv[0]))
    filename = sys.argv[1]
    property = sys.argv[2]
    main(filename, property)
