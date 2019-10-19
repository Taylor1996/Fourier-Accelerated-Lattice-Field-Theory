{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 19,
   "metadata": {},
   "outputs": [],
   "source": [
    "import re\n",
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "infile = open(\"output_hot_1M_working.txt\", \"r\")\n",
    "s = infile.read()\n",
    "tokens = s.split('[')[1:]\n",
    "\n",
    "# we discard tokens[0]\n",
    "for i,path in enumerate(tokens):\n",
    "    new_path = re.sub('\\]|\\n', '', path)\n",
    "    tokens[i] = new_path\n",
    "#x_data = [path for path in tokens[1:]]\n",
    "pathList = [[0]*len(tokens[0])]*len(tokens)\n",
    "for i in range(len(tokens)):\n",
    "    pathList[i] = [float(x_i) for x_i in tokens[i].split()]\n",
    "\n",
    "N_tau = len(pathList[0])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "120\n"
     ]
    }
   ],
   "source": [
    "# obtaining the average of the action <S>\n",
    "# to do this, we take each path and use the formula S = sum i to N_tau\n",
    "m = 1.0\n",
    "omega = 1.0\n",
    "print(N_tau)\n",
    "\n",
    "\n",
    "S_arr = np.zeros(len(pathList))\n",
    "for pNum,path in enumerate(pathList):\n",
    "    x = path\n",
    "    S = 0.0\n",
    "    for i in range(N_tau):\n",
    "        tau_plus = (i+1)%N_tau\n",
    "        tau_minus = (i-1+N_tau)%N_tau\n",
    "        #print(tau_minus)\n",
    "        x_plus = x[tau_plus]\n",
    "        x_minus = x[tau_minus]\n",
    "        S += 0.5*m*(x_plus - x[i])**2 + 0.5*m*(omega**2)*(x[i]**2)\n",
    "    S_arr[pNum] = S"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([0.20409817, 0.0854556 , 0.        , ..., 2.1135147 , 0.30194693,\n",
       "       0.55001145])"
      ]
     },
     "execution_count": 21,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# define the observable you want the error on.\n",
    "import itertools\n",
    "O = np.array(list(itertools.chain.from_iterable(pathList)))**2\n",
    "O\n",
    "# O = S_arr"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 47,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "#psi_0 = np.histogram([np.mean(path) for path in pathList],bins='auto', density=True)\n",
    "# REMEMBER THAT WE DON'T HAVE TO SQUARE psi_0 BECAUSE IN FORMING A HISTOGRAM, WE GET A PROBABILITY DENSITY FOR x, WHICH IS \n",
    "# PRECISELY WHAT |psi_0|^2 is. \n",
    "psi_0 = np.histogram([x_i for path in pathList for x_i in path],bins=100, density=True)\n",
    "# we obtain the average of each path to ploAssertionErrort psi_0 and obtain the other observables\n",
    " \n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 49,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[<matplotlib.lines.Line2D at 0x2b26e0816d8>]"
      ]
     },
     "execution_count": 49,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "# I should normalize btw\n",
    "# plt.scatter(list(psi_0[1])[1:], list(psi_0[0]**2))\n",
    "plt.scatter(list(psi_0[1])[1:], list(psi_0[0]))\n",
    "x = np.linspace(-3,3,100)\n",
    "y = (1/(np.sqrt(np.pi)))*np.exp(-x**2)\n",
    "plt.plot(x,y)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "2.718281828459045"
      ]
     },
     "execution_count": 27,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "np.exp(1)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Obtain Jacknife estimates of the errors on <x>, <x^2> \n",
    "import math\n",
    "\n",
    "\n",
    "N = len(O)\n",
    "B = int(0.01*N)\n",
    "\n",
    "############################################################################\n",
    "# if(N%B != 0.0):                                                          #\n",
    "#     # B doesn't divide N. ditch first few O_i until we have divisibility #\n",
    "#     O = O[(N/B - math.floor(N/B))* B):]                                  #\n",
    "############################################################################\n",
    "\n",
    "N_B = int(N/B)\n",
    "\n",
    "# obtain <x>, the average over the gathered configs\n",
    "tot = np.sum(O)\n",
    "avg = tot/N\n",
    "\n",
    "\n",
    "o = np.zeros(N_B)\n",
    "jack_estim = np.zeros(N_B)\n",
    "for i in range(N_B):\n",
    "    o[i] = np.sum(O[i*B:(i+1)*B])\n",
    "    jack_estim[i] = (1/(N-B)) * (tot - o[i])\n",
    "    \n",
    "#np.linalg.norm(o-avg)**2\n",
    "jack_var = ((N_B-1)/(N-B))*(np.linalg.norm(jack_estim-avg)**2)\n",
    "\n",
    "error = math.sqrt(jack_var)\n",
    "\n",
    "error"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "avg"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
