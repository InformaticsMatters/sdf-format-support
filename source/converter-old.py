#!/usr/bin/env python
import json
import re
import uuid
import datetime
from typing import List, Dict

# from rdkit import Chem

_MIME_TYPE_MAP = {'chemical/x-mdl-sdfile': 'sdf',
                  'application/x-squonk-dataset-molecule-v2+json': 'json',
                  'application/schema+json': 'json_schema'}


class ConvertFile(object):
    """Class ConvertFile

    Purpose: Converts from an input mime-type to an output mime-type.

    """

    def convert(self, from_mime_type: str, to_mime_type: str, infile: str, outfile: str) -> bool:
        """Dispatch method"""
        from_type: str = _MIME_TYPE_MAP.get(from_mime_type)
        to_type: str = _MIME_TYPE_MAP.get(to_mime_type)

        method_name = 'convert_' + str(from_type) + '_to_' + str(to_type)
        # Get the method from 'self'. Default to a lambda.
        method = getattr(self, method_name, lambda: "Invalid conversion")
        # Call the method as we return it
        return method(infile, outfile)  # pylint: disable=too-many-function-args

    # Supporting functions for Json
    def add_sdf_property(self, properties: Dict[str, str], propname: str, propvalue: List[str]):
        """"
        Adds a property (name, value) pair to the properties dictionary
        """

        if propvalue[len(propvalue) - 1]:
            # property block should end with an empty line, but some SDFs are buggy
            properties[propname] = '\n'.join(propvalue)
        else:
            properties[propname] = '\n'.join(propvalue[:-1])

    # def gen_smiles(molblock):
    #     mol = Chem.MolFromMolBlock(molblock)
    #     if mol.GetConformer(-1).Is3D():
    #         Chem.AssignAtomChiralTagsFromStructure(mol)
    #
    #     isosmiles = Chem.MolToSmiles(mol)
    #     stdsmiles =  Chem.MolToSmiles(mol, isomericSmiles=False)
    #     return isosmiles, stdsmiles

    def convert_sdf_to_json(self, infile, outfile) -> bool:
        """converts the given SDF file into a Squonk json file.
           Returns True if file successfully converted.
        """

        errors: int = 0
        lines: int = 0
        records: int = 0

        file = open(infile, 'rt')
        outfile = open(outfile, 'w')
        outfile.write('[')

        pattern: str = '^>  <(.*)>'

        try:
            molecule = {}
            molblock = []
            properties: Dict[str, str] = {}

            # Loop through molecules
            while True:
                lines += 1

                # Get next line from file
                line = file.readline()

                # if line is empty
                # end of file is reached
                if not line:
                    break

                # remove newline chars
                text = line.rstrip('\n\r')

                molblock.append(text)
                # 'M  END' signifies the end of the molblock. The properties will follow
                if text == 'M  END':
                    # print('End molblock')
                    propname = None
                    propvalue = []
                    if molblock[0]:
                        molecule['name'] = molblock[0]
                    molblock = '\n'.join(molblock)
                    molecule['molblock'] = molblock
                    # isosmiles, stdsmiles = gen_smiles(molblock)
                    # molecule['isosmiles'] = isosmiles
                    # molecule['stdsmiles'] = stdsmiles

                    # Loop through properties
                    while True:
                        lines += 1
                        line = file.readline()

                        # if line is empty
                        # end of file is reached
                        if not line:
                            break

                        text = line.strip()

                        # '$$$$' signifies the end of the record
                        if text == '$$$$':
                            # print('End record')
                            records += 1
                            if propname:
                                self.add_sdf_property(properties, propname, propvalue)

                            record = {'uuid': str(uuid.uuid4()), 'molecule': molecule, 'values': properties}
                            json_str = json.dumps(record)
                            if records > 1:
                                outfile.write(',')
                            outfile.write(json_str)

                            molecule = {}
                            molblock = []
                            properties: Dict[str, str] = {}
                            break
                        else:
                            result = re.match(pattern, text)
                            if result:
                                if propname:
                                    self.add_sdf_property(properties, propname, propvalue)
                                propname = result.group(1)
                                propvalue = []
                            else:
                                propvalue.append(text)

            outfile.write(']')

        finally:
            outfile.close()

        print(lines, errors)
        if errors > 0:
            return False
        else:
            return True

    # Supporting function for json_schema
    def is_number(self, value: str) -> int:
        """"
        Determines the number type for the property value
        :return 0 for text, 1 for float, 2 for integer
        """
        try:
            f = float(value)
        except ValueError:
            return 0
        else:
            if f.is_integer():
                return 2
            else:
                return 1

    def convert_sdf_to_json_schema(self, infile, outfile) -> bool:
        """converts the given SDF file into a Squonk json schema.
           Returns True if file successfully converted.
        """

        # sdf properties for internal json format - used in the json schema on internal conversion.
        _json_properties_sdf = {
            'uuid': {'type': 'string', 'description': 'Unique UUID'},
            'molecule': {'type': 'object', 'properties': {
                'name': {'type': 'string', 'description': 'Molecule name'},
                'molblock': {'type': 'string', 'description': 'Molfile format 2D or 3D representation'},
                'isosmiles': {'type': 'string', 'description': 'Isomeric SMILES representation'},
                'stdsmiles': {'type': 'string', 'description': 'Standardised non-isomeric SMILES representation'}
            }}}
        _lookups = {0: 'string', 1: 'number', 2: 'integer'}

        errors: int = 0
        lines: int = 0
        records: int = 0

        # Extract the properties from the file
        file = open(infile, 'rt')

        pattern = '^>  <(.*)>'
        try:
            properties = {}
            prop_types = {}
            # Loop through molecules
            while True:
                lines += 1

                # Get next line from file
                line = file.readline()

                # if line is empty
                # end of file is reached
                if not line:
                    break

                # remove newline chars
                text = line.rstrip('\n\r')
                # print(text)

                # 'M  END' signifies the end of the molblock. The properties will follow
                if text == 'M  END':
                    propname = None
                    propvalue = []

                    # Loop through properties
                    while True:
                        lines += 1
                        line = file.readline()

                        # if line is empty
                        # end of file is reached
                        if not line:
                            break

                        text = line.strip()

                        # '$$$$' signifies the end of the record
                        if text == '$$$$':
                            # print('End record')
                            records += 1
                            if propname:
                                self.add_sdf_property(properties, propname, propvalue)

                            for name in properties:
                                if name not in prop_types:
                                    prop_types[name] = 2
                                t = self.is_number(properties[name])
                                if t < prop_types[name]:
                                    prop_types[name] = t
                                # print(name, types[name])

                            properties = {}
                            break
                        else:
                            result = re.match(pattern, text)
                            if result:
                                if propname:
                                    self.add_sdf_property(properties, propname, propvalue)
                                propname = result.group(1)
                                propvalue = []
                            else:
                                propvalue.append(text)

        finally:
            file.close()

        # Add the properties from the file to the schema
        values = {}
        _json_properties_sdf['values'] = {'type': 'object', 'properties': values}
        for prop_name in prop_types:
            prop_type = prop_types[prop_name]
            values[prop_name] = {'type': _lookups[prop_type]}

        schema_sdf = {'$schema': 'http://json-schema.org/draft/2019-09/schema#',
                      'description': 'Automatically created from ' + infile + ' on ' + str(datetime.datetime.now()),
                      'properties': _json_properties_sdf}

        json_str: str = json.dumps(schema_sdf)
        with open(outfile, 'w') as schemaFile:
            schemaFile.write(json_str)

        print(lines, errors)
        if errors > 0:
            return False
        else:
            return True


if __name__ == "__main__":
    a = ConvertFile()
    print(a.convert('chemical/x-mdl-sdfile', 'application/x-squonk-dataset-molecule-v2+json', 'poses.sdf', 'out.json'))
    print (a.convert('chemical/x-mdl-sdfile', 'application/schema+json', 'poses.sdf', 'out.schema.json'))
