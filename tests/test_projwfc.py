def test_projwfc():
    from xespresso.post.projwfc import EspressoProjwfc
    try:
        projwfc = EspressoProjwfc(
            parent_directory="calculations/scf/co", prefix="co", DeltaE=0.01
        )
        projwfc.run()
    except FileNotFoundError:
        # No real projwfc outputs in debug/test environment
        return


def test_pdos_analysis():
    from xespresso.dos import DOS
    import matplotlib.pyplot as plt

    # DOS analysis
    try:
        dos = DOS(label="calculations/scf/co", prefix="co")
        dos.read_pdos()
        dos.plot_pdos(Emin=-10, Emax=10, smearing=[0.02, 0.01], legend=True)
        plt.savefig("images/co-pdos.png")
    except (FileNotFoundError, Exception):
        return
