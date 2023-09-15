"""
MicrobeRX is a tool for enzymatic reaction-based metabolite prediction in the gut microbiome.

The tool allows you to:
- Load and process reaction rules and evidences from various sources
- Generate and rank metabolite candidates for a given input compound and a set of microbes
- Visualize and explore the results using interactive plots and tables
"""

# Handle versioneer
from ._version import get_versions

# Add imports here

from .RuleGenerator import *

from .DataFiles import *

from .MetabolitePredictor import *

from .MetaboliteAnalyzer import *

from .MetaboliteVisualizer import *

__version__ = get_versions()["version"]

del get_versions

__documentation_web__ = "https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX"
__github_web__ = "https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX"
__github_issues_web__ = __github_web__ + "/issues"
