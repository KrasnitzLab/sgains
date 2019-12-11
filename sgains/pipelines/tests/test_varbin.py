import os
import numpy as np

from sgains.pipelines.varbin_pipeline import VarbinPipeline
import pytest


# @pytest.mark.skip(reason="no such file available. Need new way of testing")
@pytest.mark.parametrize('bamfile,rows', [
    # 'CJA0918.rmdup.bam',
    ('CJA0918.100.rmdup.bam', 99),
    ('CJA0918.200.rmdup.bam', 199),
])
def test_bin_counts_0918(tests_config, bamfile, rows, varbin0918):
    pipeline = VarbinPipeline(tests_config)
    bam_filename = os.path.join(
        os.environ.get("SGAINS_DATA"),
        f"CJA0918/{bamfile}"
    )

    # df = pipeline.varbin("data/test_study/CJA0918.jude58.rmdup.bam")
    df = pipeline.varbin(bam_filename)
    print(df.head())
    print(varbin0918.head())

    assert np.all(varbin0918.columns == df.columns)
    assert np.all(varbin0918.chrom == df.chrom)
    assert np.all(varbin0918.chrompos == df.chrompos)

    assert np.all(
        np.abs(varbin0918.iloc[:rows, 3] - df.iloc[:rows, 3]) <= 1
    )
    # assert np.all(
    #     np.abs(varbin0918.ratio - df.ratio) <= 1E-2
    # )
