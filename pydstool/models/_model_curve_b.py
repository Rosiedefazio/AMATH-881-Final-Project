import os

import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(0.0, 8.0 * np.pi, 700)
x = np.cos(t) * np.exp(-0.03 * t)
y = np.sin(t) * np.exp(-0.03 * t)

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(x, y, color='#d62728', linewidth=2.0)
ax.set_title('Dummy Model B')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
ax.grid(alpha=0.3)

savefig_path = os.getenv('PYDSTOOL_SAVEFIG')
if savefig_path:
    plt.savefig(savefig_path, dpi=160, bbox_inches='tight')
    print('Saved plot to %s' % savefig_path)
else:
    plt.show()
