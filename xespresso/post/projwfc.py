from xespresso.post.base import PostCalculation
from typing import Dict


class EspressoProjwfc(PostCalculation):
    """
    Interface to projwfc.x for projected density of states calculations.
    
    Projects wavefunctions onto orthogonalized atomic wavefunctions,
    calculates Lowdin charges, spilling parameter, projected DOS
    (separated into up and down components for spin-polarized systems).
    
    QE 7.5 Reference: https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html
    """

    package = "projwfc"
    package_parameters = {
        "PROJWFC": [
            "prefix",
            "outdir",
            "ngauss",
            "degauss",
            "Emin",
            "Emax",
            "DeltaE",
            "lsym",
            "diag_basis",
            "pawproj",
            "filpdos",
            "filproj",
            "filowdin",
            "lwrite_overlaps",
            "lbinary_data",
            "kresolveddos",
            "tdosinboxes",
            "n_proj_boxes",
            "irmin",
            "irmax",
            "plotboxes",
        ]
    }

    def __init__(
        self, 
        parent_directory, 
        prefix, 
        queue=False, 
        parallel="", 
        debug=False,
        dry_run=False,
        **kwargs
    ) -> None:
        """
        Initialize EspressoProjwfc.
        
        Parameters:
        - parent_directory: Directory with SCF output
        - prefix: Prefix for SCF calculation (reads prefix.save)
        - queue: Scheduler configuration
        - parallel: Parallelization flags for projwfc.x
        - debug: Enable debug logging
        - dry_run: If True, only generate input without execution
        - **kwargs: Additional projwfc parameters
        """
        super().__init__(
            parent_directory=parent_directory,
            prefix=prefix, 
            queue=queue, 
            parallel=parallel,
            debug=debug,
            dry_run=dry_run,
            **kwargs
        )
        
        # Set default values for common parameters
        if 'lsym' not in self.parameters:
            self.parameters['lsym'] = '.false.'
        if 'pawproj' not in self.parameters:
            self.parameters['pawproj'] = '.false.'
    
    def write_package_input(self):
        """Write projwfc input file with proper formatting.
        
        String parameters get quotes, logical parameters don't.
        """
        import os
        filename = os.path.join(self.directory, "%s.%si" % (self.prefix, self.package))
        defaults = self.get_defaults()
        
        with open(filename, "w") as f:
            for section, parameters in self.package_parameters.items():
                if section != "LINE":
                    f.write("&%s\n" % section)
                    for key, value in self.parameters.items():
                        if key in parameters:
                            # Skip if matches default
                            if key in defaults and defaults[key] == value:
                                continue
                            
                            # Logical parameters: .true., .false.
                            if value in ['.true.', '.false.']:
                                f.write("  {0:15s} = {1}\n".format(key, value))
                            # Numeric parameters
                            elif isinstance(value, (int, float)):
                                f.write("  {0:15s} = {1}\n".format(key, value))
                            # String parameters
                            else:
                                f.write('  {0:15s} = "{1}"\n'.format(key, value))
                    f.write("/ \n")
    
    def get_defaults(self) -> Dict:
        """Get default values for projwfc parameters per QE 7.5 documentation.
        
        Only parameters that match defaults will be skipped from the input file.
        """
        return {
            'ngauss': 0,
            'degauss': 0.0,
            'lsym': '.false.',
            'diag_basis': '.false.',
            'pawproj': '.false.',
            'lwrite_overlaps': '.false.',
            'lbinary_data': '.false.',
            'kresolveddos': '.false.',
            'tdosinboxes': '.false.',
            'n_proj_boxes': 1,
            'plotboxes': '.false.',
        }
