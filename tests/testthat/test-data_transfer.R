test_that("data transfer", {
    data.out = generate_data(n=2000, p=500, ntask=5, K=6, cl.sep=1, sig=1, batch.effect.sig=1, alpha=0.5)
    ref.data = 1:4
    int.out = CFITIntegrate(X.list=data.out$X.list[ref.data], r=6, max.niter = 10, seed=0, verbose=F)

    tf.out = CFITTransfer(Xtarget=data.out$X.list[[5]], Wref=int.out$W, max.niter = 10, seed=0, verbose=F)

    # integration, run 5 iterations for testing
    transfer.out = expect_error(CFITTransfer(Xtarget=data.out$X.list[[5]], Wref=int.out$W,
                                             max.niter = 10, seed=0, verbose=F), NA)
    transfer.out = expect_error(CFITTransfer(Xtarget=data.out$X.list[[5]], Wref=int.out$W,
                                             max.niter = 10,future.plan='multicore', seed=0, verbose=F), NA)
    transfer.out = expect_error(CFITTransfer(Xtarget=data.out$X.list[[5]], Wref=int.out$W,
                                             max.niter = 10,future.plan='multisession', seed=0, verbose=F), NA)
})

