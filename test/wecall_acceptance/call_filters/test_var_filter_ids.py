# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestVarFilterIDs(BaseTest):
    def test_should_error_with_not_allowed_var_filter_id(self):
        svc = SVCDriver(self) \
            .with_var_filters("JONNY", "SB") \
            .with_verbosity(0)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
        )
        expect = svc.call(expected_success=False)
        expect.incorrect_var_ids_error("JONNY")

    def test_should_error_with_not_allowed_var_filter_ids(self):
        svc = SVCDriver(self) \
            .with_var_filters("JONNY", "ANDY", "SB", "EDWARD", "STEFANIE") \
            .with_verbosity(0)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
        )
        expect = svc.call(expected_success=False)
        expect.incorrect_var_ids_error("JONNY", "ANDY", "EDWARD", "STEFANIE")
