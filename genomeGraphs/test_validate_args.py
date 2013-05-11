#!/usr/bin/env py.test

import validate_args as va

def test_validate_ymax():
    ymax= ['max']
    assert va.validate_ymax(ymax) is True
    ymax= ['indiv', 1, '10.2']
    assert va.validate_ymax(ymax) is True
    ymax= ['max', 1]
    assert va.validate_ymax(ymax) is False
    ymax= ['max', 'indiv']
    assert va.validate_ymax(ymax) is False
    ymax= ['indiv', 1, 'foo']
    assert va.validate_ymax(ymax) is False
