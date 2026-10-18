"""xs_locate.py with the USGS Cascadia CVM v1.7 regional models registered (USGS_OR etc.; see lfe_array_discovery/scripts/ab_hypo_cvm.py)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import xs_locate as L        # puts lfe_array_discovery/scripts on sys.path and imports ab_hypo as H
import ab_hypo_cvm           # noqa: F401  registers USGS_* models into ab_hypo.MODELS
if __name__ == "__main__":
    L.main()
