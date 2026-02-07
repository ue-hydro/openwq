Join us!
========

Would you like to contribute to OpenWQ?

Do you need help to couple OpenWQ to your model?

Do you have an interesting idea involving OpenWQ?

Contact us through our GitHub page (``ue-hydro``) at the University of Evora: https://github.com/ue-hydro


How to contribute
~~~~~~~~~~~~~~~~~~

OpenWQ welcomes contributions from the community. Here are some ways to get involved:

**Report bugs**: Open an issue on the `GitHub repository <https://github.com/ue-hydro>`_ describing the problem, including the configuration files and error messages.

**Suggest features**: Open a GitHub issue describing the feature you would like to see, including use cases and expected behavior.

**Submit code**: Fork the repository, create a feature branch, and submit a pull request. Please ensure your code:

* Compiles without warnings in both ``debug`` and ``fast`` modes
* Passes the existing :doc:`synthetic tests <synthetic_tests>`
* Includes appropriate comments describing the purpose of new functions
* Follows the existing code style (indentation, naming conventions)

**Add a new host model coupling**: Follow the :doc:`Coupler Guide <5_3_Coupler_guide>` to integrate OpenWQ with your hydrological model. The process requires implementing four API calls in your model's source code.

**Improve documentation**: This documentation is written in reStructuredText and built with Sphinx. Contributions to improve clarity, fix errors, or add examples are always welcome.


Development setup
~~~~~~~~~~~~~~~~~~

1. Clone the repository::

    git clone --recurse-submodules https://github.com/ue-hydro/openwq.git

2. Build in debug mode for development::

    cmake -DHOST_MODEL_TARGET=openwq -DCMAKE_BUILD_TYPE=debug .
    make -j 4

3. Run the synthetic tests to verify your setup::

    # Navigate to the synthetic test directory
    # Run each test case and compare against analytical solutions

4. Make your changes on a feature branch and submit a pull request.
