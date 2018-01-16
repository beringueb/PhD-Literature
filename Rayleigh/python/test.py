import numpy as np
import tensorflow as tf

FLAGS = NONE

def main(_):
   #Paramters
   l_max=500
   iteration=10000
   learning_rate=1

   #Files loc
   dir='/home/bb510/Data/Test/clsfile_0256_0500_'
   filename=[(dir+'%04d.txt' %i) for i in range(10000)]
   
   #Variables and constants
   PI=tf.Variable(tf.zeros([501,501],dtype=tf.float32),dtype=tf.float32)
x=tf.placeholder(tf.float32, [l_max+1,1])
y=tf.matmul(PI,x)
y_true=tf.placeholder(tf.float32,[l_max+1,1])
cross_entropy=tf.reduce_sum(tf.nn.softmax_cross_entropy_with_logits(labels=y_true, logits=y))
train_step=tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)
sess=tf.InteractiveSession()
init=tf.initialize_all_variables
for i in range(500) :
   l,cl_out,cl_in,cl_true=np.loadtxt(filename[i], usecols=(0,1,2,3),unpack=True,dtype=float)
   sess.run(train_step,feed_dict={x:np.float32(cl_in.reshape(l_max+1,1)),y_true:np.float32(cl_out.reshape(l_max+1,1))})

