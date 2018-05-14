Simple help functions for testing bash scripts
==============================================

Download:

```
wget https://github.com/dariober/bioinformatics-cafe/blob/master/bashTestFunctions/bashTestFunctions.sh?raw=true -O bashTestFunctions.sh
```

Example:

```
#!/bin/bash

source bashTestFunctions.sh

pprint 'Test exit code'
echo "FOOBAR"
assertEquals 0 $?
```

See script itself for available functions.
