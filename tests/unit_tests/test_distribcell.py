import pytest

import openmc


def test_invalid_cell():
    pwr_pincell = openmc.examples.pwr_pin_cell()


    # add a distribcell tally with a cell that doesn't exist in the model
    cell = openmc.Cell()

    dcell_filter = openmc.DistribcellFilter(cell)

    tally = openmc.Tally()
    tally.filters = [dcell_filter]
    tally.scores = ['flux']

    pwr_pincell.tallies = [tally]

    with pytest.raises(RuntimeError, match=r'Could not find cell'):
        pwr_pincell.run()


if __name__ == '__main__':
    pytest.main()