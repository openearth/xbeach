Sphinx cheatsheet
=================

Sphinx source files are plain text files.  You can simply write your
text in a text editor of your choice. Separate paragraphs by a blank
line and you will do ok.

Some simple markup code is available to set text formatting, insert
headers, figures and tables and make cross-references throughout the
document. The most important markup code is listed below.

Text formatting
---------------

.. code-block:: rst
                
   Simple text formatting is available:

   *italic*
   **bold**

   Bullet lists:
   * bullet 1
   * bullet 2
   * bullet 3

   Numbered lists:
   #. item 1
   #. item 2
   #. item 3

Headers
-------

.. code-block:: rst

   Headers can be added by using different underlining. Optionally, an
   header can have a label. You can make references to a header using
   its label as follows: Section :ref:`sec-heading`.

   Heading 1
   =========

   Heading 2
   ---------

   Heading 3
   ^^^^^^^^^

   Heading 4
   ~~~~~~~~~

   .. _sec-heading:

   Heading with label
   ------------------

Figures
-------

.. code-block:: rst

   Figure are created using the following code. A figure points to a
   filename and has a caption. Optionally, you can add a label and set
   the width, alignment and other common figure properties. You can
   make references to a figure using its label as follows: Figure
   :numref:`fig-coordsys-rect`.

   .. _fig-coordsys-rect:

   .. figure:: image8.png
      :width: 400px
      :align: center

      Rectangular coordinate system of XBeach

Tables
------

.. code-block:: rst

   Tables are created using the following code. A table has a
   captionfigure points to a filename and has a caption. Optionally,
   you can add a label. You can make references to a table using its
   label as follows: Table :numref:`tab-wave-breaking`.

   .. _tab-wave-breaking:

   .. table:: Different wave breaking formulations implemented

      +-----------------------------+-----------------+------------------+
      | Wave breaking formula       | Type of waves   | keyword          |
      +=============================+=================+==================+
      | Roelvink (1993a)            | Instationary    | roelvink1        |
      +-----------------------------+-----------------+------------------+
      | Roelvink (1993a) extended   | Instationary    | roelvink2        |
      +-----------------------------+-----------------+------------------+
      | Daly et al. (2010)          | Instationary    | roelvink_daly    |
      +-----------------------------+-----------------+------------------+
      | Baldock et al. (1998)       | Stationary      | baldock          |
      +-----------------------------+-----------------+------------------+
      | Janssen & Battjes (2007)    | Stationary      | janssen          |
      +-----------------------------+-----------------+------------------+
      
Equations
---------

.. code-block:: rst

   Equations can be created using Latex math code or by copying and
   pasting from MathType. An equation is created using the following
   code. An equation has a label. You can make references to an
   equation using its label as follows: Equation :eq:`eq-eikonal`.

   .. math::
      :label: eikonal

      \begin{array}{l}
      {k_{x} =k_{x}^{n-1} +k_{x}^{:} } \\
      {k_{y} =k_{y}^{n-1} +k_{y}^{:} } \\
      \end{array}

   Inline math can be added as follows: :math:`k_x`.

Citations
---------

.. code-block:: rst

   Literature references are added using a BibTeX file maintained
   through the Mendeley XBeach group. You can cite a reference using
   :cite:`Roelvink2009`.

Parameter lists
---------------

.. code-block:: rst

   Parameter lists can be exported using the XBeach maintainance
   toolbox. These exports have syntax according to the following
   code. This code is rendered as a normal table. References to
   individual parameters can be made as follows: :par:`<Brfac>`.

   .. partable:: Overview of available keywords

      ARC
        :advanced:
        :description:   Switch for active reflection compensation at seaward boundary
        :units:         -
        :default:       1
        :range:         0 - 1
      Brfac
        :advanced:
        :description:   Calibration factor surface slope
        :units:         -
        :default:       1.0
        :range:         0.0 - 1.0
      C
        :description:   Chezy coefficient
        :units:         m^0.5s^-1
        :default:       55.0
        :range:         20.0 - 100.0

Source code references
----------------------

.. code-block:: rst

   References to the Fortran source code can be made with the
   following code. Note that the references are made to the actual
   name of the module, which is not necessarily the same as the
   filename.

   Module: :f:mod:`flow_timestep_module`

   Function: :f:mod:`flow_timestep_module/visc_smagorinsky`
      
Includes
--------

.. code-block:: rst

   Other RST files can be included with the following code. Including
   files is used to embed the automatically generated parameter lists
   into the documentation.

   .. include:: tables/partable_all.tab

File contents
-------------

.. code-block:: rst

   Literal file contents can be added using the following code.

   **jonswap.txt**

   .. code-block:: text

      Hm0 = 0.8
      Tp = 8
      mainang = 285.
      gammajsp = 3.3
      s = 10.
      fnyq = 0.3

Footnotes and URL's
-------------------

.. code-block:: rst

   References to footnotes and URL's are made using the following code:

   Footnote: [#1]_

   URL: [http://www.xbeach.org]_

   The footnote itself is made as follows:

   .. [#1] This is a footnote
           
See also
--------

.. code-block:: rst

   Yellow frames with references to additional information, for
   example Fortran source code descriptions, can be made using the
   following code.

   .. seealso:: Wave shapes are implemented in :f:func:`morphevolution/RvR`
                and :f:func:`morphevolution/vT`.
                
