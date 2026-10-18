import numpy as np
from codex_lfe.network_probability_scan import score_grid


def test_travel_alignment_and_nonwrapping_shift_sign():
    # Four arrivals from a known origin must align at that origin. Positive
    # synthetic time shifts move arrivals later, rather than wrapping the day.
    offsets = np.array([[3, 5, 8, 11], [4, 4, 7, 12]])
    scores = np.zeros((4, 160), np.float32)
    for station, delay in enumerate(offsets[0]):
        scores[station, 70 + delay] = 5
    indices = np.arange(20, 120)
    all_scores, gated, node = score_grid(scores, offsets, indices, np.zeros(4, int), 'cpu')
    assert indices[np.argmax(gated)] == 70
    assert node[np.argmax(gated)] == 0
    assert gated.max() == 5
    assert all_scores.max() == 5
    _, shifted, _ = score_grid(scores, offsets, indices, np.full(4, 9), 'cpu')
    assert indices[np.argmax(shifted)] == 79
    _, broken, _ = score_grid(scores, offsets, indices, np.array([0, 9, 18, 27]), 'cpu')
    assert np.isneginf(broken).all()
