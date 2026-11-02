import sitoolz.primebase as pb


def test_factorize_and_largest():
    assert pb.factorize(360) == {2: 3, 3: 2, 5: 1}
    assert pb.factorize(259313) == {257: 1, 1009: 1}
    assert pb.largest_factor(259313) == 1009
    assert pb.largest_factor(1009732533765288) == 3497381   # N-1 is smooth


def test_pi_exact():
    assert pb.pi(10) == 4
    assert pb.pi(100) == 25
    assert pb.pi(1_000_000) == 78_498
    assert pb.pi(10 ** 8) == 5_761_455       # from the stored base (searchsorted)


def test_pi_beyond_base_via_lucy():
    # 10^9 is past the stored primes -> exact Lucy count
    assert pb.pi(10 ** 9) == 50_847_534


def test_count_range_and_octaves():
    assert pb.count_range(2, 100) == 25
    assert pb.count_range(10 ** 12, 10 ** 12 + 100) == 4
    oc = pb.octaves(10, 12)
    assert oc[0]["exact"] == 137            # primes in [1024, 2048)
    assert oc[0]["li_est"] > 130


def test_li_close_to_pi():
    assert abs(pb.li(10 ** 6) - pb.pi(10 ** 6)) < 200


def test_neighborhood_anchor():
    info = pb.neighborhood(1009732533765289, width=6)
    assert info["prime"] is True
    # N+3 = 4 * a giant prime -> a giant offset exists as a proof anchor
    assert info["anchor"] is not None
    assert info["anchor"]["giant"] is True


def test_rollball_descend_and_up():
    chain = pb.descend(9_999_999_999_999_917)
    assert chain[0] == 9_999_999_999_999_917
    assert chain[-1] == 73                   # rolls down to seed 73
    ups = pb.roll_up(73, count=3)
    assert all(pb.is_prime(N) for _, N, _ in ups)


def test_descend_refuses_past_reach():
    import pytest
    with pytest.raises(ValueError):
        pb.descend(pb.REACH * 2)             # stays pure, no guessing


def test_pocklington_prove_and_certificate():
    # matches L.is_pocklington_prime(N, a): single-witness call, truthy on prime
    assert pb.is_pocklington_prime(29, 7)
    assert pb.is_pocklington_prime(2 * 29 + 1, 29)
    # shape-but-not-prime: 13817611737 = 2^3*1727201467+1, composite
    assert pb.is_pocklington_prime(13817611737) is False
    assert pb.is_pocklington_prime(14488816045513223)["q"]  # a real chain prime
    # full certificate for a 31-digit prime, verified from scratch
    P = 1019559789743953361505824483047
    cert = pb.certificate(P)
    assert isinstance(cert, list) and pb.verify_certificate(cert)
    assert cert[0]["q"] == 10698977642264269


def test_self_prove():
    # cascade: descend to a base prime, climb the chain to N, witness from below
    N = 10195597897439533615058244830471234324417
    chain = pb.self_prove(N)
    assert isinstance(chain, list)
    assert chain[-1]["N"] == N                       # top of the ladder is N
    assert chain[0].get("proof") == "base"           # rooted in a base-proven prime
    # Every rung is either proven outright by the base or witnessed by the rung below
    # it. Which rungs fall in which camp depends on REACH -- a bigger base proves more
    # of the ladder directly -- so assert the structure, not a particular witness.
    for lower, upper in zip(chain, chain[1:]):
        assert upper.get("proof") == "base" or upper["a"] == lower["N"]
    assert pb.self_prove(1009732533765251 * 1009732533765211) is False


def test_roll_chain_returns_all_trues():
    assert pb.roll_chain(2) == [2, 5, 11, 23, 47]
    assert pb.roll_chain(89) == [89, 179, 359, 719, 1439, 2879]
    assert all(pb.is_prime(x) for x in pb.roll_chain(2))


def test_chain_descend_verifies_rungs():
    links = pb.chain_descend(47)             # Cunningham chain 2-5-11-23-47
    assert [L["x"] for L in links] == [47, 23, 11, 5]
    assert links[0]["rung"] and links[1]["rung"] and links[2]["rung"]
    # a composite's odd parts aren't prime -> no genuine rungs
    N = 1009732533765211 * 1009732533765251
    assert all(not L["rung"] for L in pb.chain_descend(N)[:5])
