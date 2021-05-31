# Tests

The following tests can be run using: ./test/runtests.sh

## success/1

### Test 1.1: SDF import - UUID added

Generate file where UUID is added. 

### Test 1.2 1.2: SDF import - No UUID 

Generate file with UUID not overidden. 
This will produce error if the input file has mol blocks that do not already have a UUID in the first line.  

### Test 1.3: SDF import - UUID added and Name rewritten

Generate file with UUID and any existing molecule name written to a new parameter.

### Test 1.4: SDF import - Invalid UUID in header

Generate file where there is a record with an invalid UUID.

### Test 1.5: SDF convert to json

Here we covert our test/success/**1** dataset to a JSON file...

### Test 1.6: SDF convert to json schema

Here we covert our test/success/**1** dataset to a JSON schema file...

## failure/1

### Test 1.1: SDF import - UUID added ad Name not rewritten

Generate file no UUID and any existing molecule name written to a new parameter - invalid combination.

