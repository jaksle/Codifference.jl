# Codifference

[![Build Status](https://github.com/jaksle/Codifference.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jaksle/Codifference.jl/actions/workflows/CI.yml?query=branch%3Amaster)

These package contains tools for estimating and analysing dispersion and dependence which can be calculated for both finite and infinite variance models. The main implemented functions are:
- Log characteristic function
  
$$l^\theta(X) = -\frac{2}{\theta^2}\ln\mathbb{E}\mathrm{e}^{\mathrm{i}\theta X}$$

available as `lcf(X, θ)`.
- Symmetric codifference
  
$$s^\theta(X,Y) = \frac{1}{4}\left(l^\theta(X+Y)-l^\theta(X-Y)\right)=\frac{1}{2\theta^2}\ln\frac{\mathbb{E}\mathrm{e}^{\mathrm{i}\theta(X-Y)}}{\mathbb{E}\mathrm{e}^{{\mathrm{i}\theta(X+Y)}}}$$

available as `cod(X, Y, θ, :s)`.
- Two types of asymmetric codifference

$$c_\pm^\theta (X,Y) = \pm\frac{1}{2}\left(l^\theta(X)+l^\theta(Y)-l^\theta(X\mp Y)\right) =\pm\frac{1}{\theta^2}\ln\frac{\mathbb{E}\mathrm{e}^{\mathrm{i}\theta(X\mp Y)}}{\mathbb{E}\mathrm{e}^{\mathrm{i}\theta X}\mathbb{E}\mathrm{e}^{\mathrm{i}\theta Y}}$$

available as `cod(X, Y, θ, :+)` or `cod(X, Y, θ, :-)`.

Moreover, the large sample asymptotic distributions of these quantities are available using `lcfAsymptDist` and `codAsymptDistr`. Asymptotic confidence invervals can be calculated using `lcfConfInterval` and `codConfInterval`.
