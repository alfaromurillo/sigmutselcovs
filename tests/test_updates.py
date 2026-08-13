"""Update-check tests (mocked HTTP)."""

import json


import sigmutselcovs.updates as updates
from sigmutselcovs.updates import check_updates, location_sources


def test_sources_json_parses_and_methods_exist():
    raw = json.loads(location_sources.read_text())
    assert raw["schema_version"] == 1
    for name, entry in raw["sources"].items():
        assert entry["check"]["method"] in updates._METHODS, name


class _HeadResponse:
    def __init__(self, headers):
        self.headers = headers

    def raise_for_status(self):
        pass


class _HeadSession:
    def __init__(self, headers, fail=False):
        self.headers = headers
        self.fail = fail

    def head(self, url, timeout=None, allow_redirects=True):
        if self.fail:
            raise OSError("connection refused")
        return _HeadResponse(self.headers)


def _one_source(tmp_path, known):
    raw = {
        "schema_version": 1,
        "checked": None,
        "sources": {
            "gtex": {
                "description": "d",
                "url": "http://x/f.gct.gz",
                "check": {
                    "method": "http_head",
                    "compare": ["content-length"],
                },
                "known": known,
            }
        },
    }
    path = tmp_path / "sources.json"
    path.write_text(json.dumps(raw))
    return path


def test_status_ok_changed_unknown(tmp_path):
    session = _HeadSession({"content-length": "123"})

    path = _one_source(tmp_path, {"content-length": "123"})
    frame = check_updates(sources_path=path, session=session)
    assert frame.loc[0, "status"] == "ok"

    path = _one_source(tmp_path, {"content-length": "999"})
    frame = check_updates(sources_path=path, session=session)
    assert frame.loc[0, "status"] == "changed"

    path = _one_source(tmp_path, {})
    frame = check_updates(sources_path=path, session=session)
    assert frame.loc[0, "status"] == "unknown"


def test_unreachable_does_not_raise(tmp_path):
    path = _one_source(tmp_path, {"content-length": "123"})
    frame = check_updates(
        sources_path=path, session=_HeadSession({}, fail=True)
    )
    assert frame.loc[0, "status"] == "unreachable"
    assert "connection refused" in frame.loc[0, "notes"]


def test_update_file_rewrites_known(tmp_path):
    path = _one_source(tmp_path, {})
    session = _HeadSession({"content-length": "123"})
    check_updates(
        sources_path=path, session=session, update_file=True
    )
    raw = json.loads(path.read_text())
    assert raw["sources"]["gtex"]["known"] == {
        "content-length": "123"
    }
    assert raw["checked"] is not None
    # next run is ok
    frame = check_updates(sources_path=path, session=session)
    assert frame.loc[0, "status"] == "ok"


def test_sources_filter(tmp_path):
    path = _one_source(tmp_path, {})
    frame = check_updates(
        sources_path=path, sources=["nope"], session=_HeadSession({})
    )
    assert frame.empty
