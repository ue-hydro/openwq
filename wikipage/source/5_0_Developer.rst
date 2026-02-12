Join Us
=======

OpenWQ is developed by the `UE-Hydro group <https://github.com/ue-hydro>`_ at the University of Évora. We welcome contributions, collaborations, and new ideas!

**Get in touch:** https://github.com/ue-hydro

Ways to Contribute
~~~~~~~~~~~~~~~~~~

* **Report bugs**: Open an issue describing the problem with configuration files and error messages
* **Suggest features**: Describe your use case and expected behavior in a GitHub issue
* **Submit code**: Fork the repository, create a feature branch, and submit a pull request
* **Add host model coupling**: Follow the :doc:`Coupler Guide <5_3_Coupler_guide>` to integrate OpenWQ with your model
* **Improve documentation**: Fix errors, add examples, or improve clarity

Code Requirements
~~~~~~~~~~~~~~~~~

Pull requests should:

* Compile without warnings in ``debug`` and ``fast`` modes
* Pass the existing :doc:`synthetic tests <synthetic_tests>`
* Include comments for new functions
* Follow existing code style

Development Setup
~~~~~~~~~~~~~~~~~

::

    # Clone repository
    git clone --recurse-submodules https://github.com/ue-hydro/openwq.git

    # Build in debug mode
    cmake -DHOST_MODEL_TARGET=openwq -DCMAKE_BUILD_TYPE=debug .
    make -j 4

    # Run synthetic tests to verify setup
    # Make changes on a feature branch and submit a pull request
