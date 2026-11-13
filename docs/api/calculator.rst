ElasticConstants (``elastic_ase.calculator``)
=============================================

.. autoclass:: elastic_ase.calculator.ElasticASE
   :members: calculate, runs, status, missing, is_complete, run_pending, write_inputs, run_jobs,
             read_results, vasp_input, vasp_status, vasp_output, qe_input, qe_status, qe_output,
             plot_fits, save, load, compliance, directional_analysis, Young, shear, poisson,
             linear_compressibility, visualize, to_elastic_vis, export_matrix

   .. automethod:: __init__

The DFT codes known to ``code=`` are described in :mod:`elastic_ase.codes`.

.. autofunction:: elastic_ase.calculator.auto_detect_symmetry
.. autofunction:: elastic_ase.calculator.validate_orientation
.. autofunction:: elastic_ase.calculator.standardize_and_reorient
