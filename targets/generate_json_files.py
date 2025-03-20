import os
import json

if __name__ == "__main__":
    targets = ['ABAHIW', 'ABAKIZ', 'ABADOX', 
               'ABABIP', 'GASQOK', 'ABEKIE',
               'NIWPUE01','ABEKIF', 'APUFEX', 
               'ABEHAU', 'TITTUO', 'EGEYOG',
               'ABOBUP', 'XIDTOW', 'ACNCOB10',
               'TACXUQ', 'ACAZFE', 'NIVHEJ',
               'ADUPAS', 'DAJLAC', 'OFOWIS',
               'CATSUL', 'HESMUQ01', 'GUDQOL',
               'ABEVAG', 'AKOQOH', 'ADARUT',
               'AFECIA', 'ACOVUL', 'AFIXEV']
    
    output_dir = './targets/'
    os.makedirs(output_dir, exist_ok=True)
    
    for index, target_entry in enumerate(targets, start=1):
        print(f"Processing target: {target_entry}")
        
        json_data = {
            "task": f"{index}_{target_entry}",
            "output_dir": "./",
            "load_from_previous": False,
            "logging": False,
            "monitor_app": False,
            "termination_exit": False,
            "scoring_functions": [
                {
                    "name": "HSR",
                    "run": True,
                    "parameters": {
                        "prefix": target_entry,
                        "ref_molecule": f"configs/3D_Benchmark/targets/{index}_{target_entry}.txt",
                        "generator": "obabel",
                        "n_jobs": 8,
                        "timeout": 10,
                        "save_files": False
                    }
                }
            ],
            "scoring": {
                "metrics": [
                    {
                        "name": f"{target_entry}_HSR_score",
                        "filter": False,
                        "weight": 1.0,
                        "modifier": "raw",
                        "parameters": {}
                    }
                ],
                "method": "single"
            },
            "diversity_filter": {
                "run": False
            }
        }
        
        filename = f"{output_dir}{index}_{target_entry}.json"
        with open(filename, 'w') as f:
            json.dump(json_data, f, indent=4)
        
        print(f"Target {target_entry} processed and saved as {filename}.")
