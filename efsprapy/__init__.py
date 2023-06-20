from os import path

# make root directory of this app which will be used 1. when running the app; 2. pyinstaller at compiling the app.
if path.exists(path.dirname(__file__)):
    # this path should be used when running the app as a Python package (non compiled) and/or pyinstaller at compiling
    # stage.
    __root_dir__ = path.realpath(path.dirname(__file__))
elif path.exists(path.dirname(path.dirname(__file__))):
    # the path will become invalid when the app run after compiled as the dirname `fsetoolsGUI` will disappear.
    # instead, the parent folder of the project dir will be used.
    __root_dir__ = path.realpath(path.dirname(path.dirname(__file__)))
else:
    raise IsADirectoryError(
        f'Project root directory undefined: '
        f'{path.dirname(__file__)} nor '
        f'{path.dirname(path.dirname(__file__))}'
    )

"""
VERSION IDENTIFICATION RULES DOCUMENTED IN PEP 440 ARE FOLLOWED.

Version scheme
==============

Distributions are identified by a public version identifier which supports all defined version comparison operations

The version scheme is used both to describe the distribution version provided by a particular distribution archive, and
to place constraints on the version of dependencies needed in order to build or run the software.

Public version identifiers
--------------------------

The canonical public version identifiers MUST comply with the following scheme:

`[N!]N(.N)*[{a|b|rc}N][.postN][.devN]`

Public version identifiers MUST NOT include leading or trailing whitespace.

Public version identifiers MUST be unique within a given distribution.

See also Appendix B : Parsing version strings with regular expressions which provides a regular expression to check
strict conformance with the canonical format, as well as a more permissive regular expression accepting inputs that may
require subsequent normalization.

Public version identifiers are separated into up to five segments:

    - Epoch segment: N!
    - Release segment: N(.N)*
    - Pre-release segment: {a|b|rc}N
    - Post-release segment: .postN
    - Development release segment: .devN

"""

__version__ = "0.0.1"
