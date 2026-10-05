"""
_optforce_core - Vendored OptForce discovery algorithm (FVA -> MUST -> FORCE).

These modules (``fva``, ``must_sets``, ``optforce``, and the ``utils`` helper
``fva`` needs) are the OptForce candidate-discovery math, vendored here so the
design engine is self-contained and the standalone ``optforce`` package can be
removed without affecting it. Only the OptForce *ranking pipeline* was dropped;
this is the underlying algorithm, unchanged.
"""
