package command

import (
	"os"
	"os/exec"
)

func detectJobManager() string {
	var path string
	var err error
	if path, err = exec.LookPath("sbatch"); err == nil {
		// SLURM
		return path
	}
	if os.Getenv("SGE_ROOT") != "" || os.Getenv("GRID_HOME") != "" {
		// Grid Engine

	}
	if path, err = exec.LookPath("qsub"); err == nil {
		// TORQUE
		return path
	}
	// None / Unsupported
	return ""
}
