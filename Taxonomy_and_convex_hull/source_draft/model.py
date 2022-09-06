"""

@Author: Niksa Praljak
"""

import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Conv2DTranspose, Layer, Lambda, Conv2D, MaxPooling2D, Dense, Dropout, Input, Flatten, BatchNormalization, Activation, Reshape

from keras import backend as K
from keras.constraints import unit_norm, max_norm
from keras.layers.advanced_activations import LeakyReLU



'Original model architecture (made during the autumn quarter):'





# create MMD loss

def compute_kernel(x, y):
    """
    

    args:
        
    """
    x_size = K.shape(x)[0]
    y_size = K.shape(y)[0]
    dim = K.shape(x)[1]
    tiled_x = K.tile(K.reshape(x, [x_size, 1, dim]), [1, y_size, 1])
    tiled_y = K.tile(K.reshape(y, [1, y_size, dim]), [x_size, 1, 1])
    return K.exp(-K.mean(K.square(tiled_x - tiled_y), axis = 2) / K.cast(dim, 'float32'))

def compute_mmd(x, y):
    """
    

    args:
        x --> samples from model embedding distribution
        y --> noisy samples from a normal distribution

    return:
        MMD divergence measurement/similarity value
    """ 
    x_kernel = compute_kernel(x,x)
    y_kernel = compute_kernel(y,y)
    xy_kernel = compute_kernel(x,y)
    
    return K.mean(x_kernel) + K.mean(y_kernel) - 2 * K.mean(xy_kernel)


# create loss layer:

def custom_VAEloss(train_z, train_xr, train_x):
    """
    
    args:
        train_z --> latent space embeddings
        train_xr --> reconstructed input samples
        train_x --> true input samples

    return:
        loss_REC --> reconstruction loss
        loss_MMD --> divergence/regularization loss

    """
    # sampling from the random noise ...
    batch_size = K.shape(train_z)[0]
    latent_dim = K.int_shape(train_z)[1]
    true_samples = K.random_normal(shape = (batch_size, latent_dim), mean = 0., stddev = 1.)
    
    # calc MMD loss
    loss_MMD = compute_mmd(true_samples, train_z)
    
    # calc the reconstruction loss (i.e. negative log-likelihood)
    loss_REC = K.mean(K.square(train_xr - train_x))

    return loss_REC + 2*loss_MMD




# create model architecture:

def orig_MMD_VAE(seq_len, aa_var, intermediate_dim, latent_dim, alpha = 0.1):
    """
    
    args:
        seq_len --> protein sequence length
        aa_var --> amino acid vocabulary length
        intermediate_dim --> width of the hidden layers
        latent_dim --> size of the latent space
        alpha --> leaky relu hyperparameter

    return:
        MMD_VAE --> the whole MMD-VAE model architecture
        encoder --> the bottom component of the autoencoder
        decoder --> the top component of the autoencoder
    """


    # encoder component:
    encoder_input = Input(shape = (seq_len, aa_var))
    x = Flatten()(encoder_input)

    x = Dense(intermediate_dim*1.5, kernel_initializer = 'random_normal',)(x)
    x = LeakyReLU(alpha)(x)
    x = Dropout(0.3)(x)
    x = Dense(intermediate_dim*1.5, kernel_initializer = 'random_normal',)(x)
    x = LeakyReLU(alpha)(x)
    x = Dense(intermediate_dim*1.5, kernel_initializer = 'random_normal',)(x)
    x = LeakyReLU(alpha)(x)
    
    encoder_output = Dense(latent_dim)(x)
    
    encoder = Model(encoder_input, encoder_output, name = 'encoder')

    # decoder component:
    decoder_input = Input(shape = (latent_dim,), name = 'z_sampling')

    x = Dense(intermediate_dim*1.5, kernel_initializer = 'random_normal', )(decoder_input)
    x = LeakyReLU(alpha)(x)
    x = Dense(intermediate_dim*1.5, kernel_initializer = 'random_normal',)(x)
    x = LeakyReLU(alpha)(x)
    x = Dropout(0.7)(x)
    x = Dense(intermediate_dim*1.5, kernel_initializer = 'random_normal',)(x)
    x = LeakyReLU(alpha)(x)

    outputs = Dense(seq_len*aa_var, activation = 'linear', kernel_initializer = 'random_normal')(x)
    outputs = Reshape((seq_len, aa_var))(outputs)
    softmax = tf.keras.activations.softmax(outputs, axis = -1)

    decoder = Model(decoder_input, softmax, name = 'decoder')

    # connect the pieces
    train_z = encoder(encoder_input)
    train_xr = decoder(train_z)
    MMD_VAE = Model(encoder_input, train_xr)

    # compile the whole model
    loss = custom_VAEloss(train_z, train_xr, encoder_input)
    MMD_VAE.add_loss(loss)

    
    return MMD_VAE, encoder, decoder


