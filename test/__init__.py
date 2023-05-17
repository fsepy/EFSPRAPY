from efsprapy import __version__


def test_version_canonical(version=__version__):
    import re
    check = re.match(
        r"^([1-9][0-9]*!)?"
        r"(0|[1-9][0-9]*)"
        r"(\.(0|[1-9][0-9]*))*"
        r"((a|b|rc)(0|[1-9][0-9]*))?"
        r"(\.post(0|[1-9][0-9]*))?"
        r"(\.dev(0|[1-9][0-9]*))?$",
        version,
    )
    assert check is not None
