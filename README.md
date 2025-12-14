# Codifference

[![Build Status](https://github.com/jaksle/Codifference.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jaksle/Codifference.jl/actions/workflows/CI.yml?query=branch%3Amaster)

- Log characteristic function
  
$$l^\theta(X) = -\frac{2}{\theta^2}\ln\mathbb{E}\mathrm{e}^{\mathrm{i}\theta X}$$

available as `lcf(X, θ)`.
- Symmetric codifference
  
$$s^\theta(X,Y) = \frac{1}{4}\left(l^\theta(X+Y)-l^\theta(X-Y)\right)=\frac{1}{2\theta^2}\ln\frac{\mathbb{E}\mathrm{e}^{\mathrm{i}\theta(X-Y)}}{\mathbb{E}\mathrm{e}^{{\mathrm{i}\theta(X+Y)}}}$$

available as `cdf(X, Y, θ, :s)`.
- Two types of asymmetric codifference

$$c_\pm^\theta (X,Y) = \pm\frac{1}{2}\left(l^\theta(X)+l^\theta(Y)-l^\theta(X\mp Y)\right) =\pm\frac{1}{\theta^2}\ln\frac{\mathbb{E}\mathrm{e}^{\mathrm{i}\theta(X\mp Y)}}{\mathbb{E}\mathrm{e}^{\mathrm{i}\theta X}\mathbb{E}\mathrm{e}^{\mathrm{i}\theta Y}}$$

available as `cdf(X, Y, θ, :+)` or `cdf(X, Y, θ, :-)`.
