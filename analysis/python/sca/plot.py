import matplotlib.pyplot as plt
import numpy as np


def latent_activity(axs, t, Z, save_path=None):

    n_time, n_sen, n_comp = Z.shape

    trial_colors = plt.cm.rainbow(np.linspace(0, 1, n_sen))

    for i, ax in enumerate(axs.ravel()):
        if i >= n_comp:
            ax.axis("off")
            continue
        
        for j in range(n_sen):
            z = Z[:, j, i]
            z = -z if np.mean(z) < 0 else z
            ax.plot(t, z, color=trial_colors[j])
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(f"Component {i+1}")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Activation")
    
    if save_path is not None:
        plt.savefig(save_path)
    else:
        plt.show()