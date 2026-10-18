import pytest
from codex_lfe.map_roster import read_map_stations


def test_reads_station_metadata_without_evaluating_large_payload(tmp_path):
    path=tmp_path/'map.html'
    path.write_text('<script>const DATA={"PGC":{"lat":48.65,"lon":-123.45,"nfam":20}},DAILY=not_json_or_code_to_execute;</script>')
    assert read_map_stations(path)=={'PGC':{'latitude':48.65,'longitude':-123.45}}


def test_missing_map_payload_rejected(tmp_path):
    path=tmp_path/'map.html';path.write_text('<html>No embedded station metadata</html>')
    with pytest.raises(ValueError,match='missing'):read_map_stations(path)
