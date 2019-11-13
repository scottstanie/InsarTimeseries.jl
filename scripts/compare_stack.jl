function compare_stack(stack, nums; vm=3, dvm=2, cmap="seismic_wide_y")

    fig, axes = plt.subplots(2, length(nums), sharex=true, sharey=true)
    ref = stack[:, :, end-nums[1]]

    for (idx, n) in enumerate(nums)
        axes[1, idx].imshow(stack[:,:,end - n], vmax=vm, vmin=-vm, cmap=cmap)
        axes[2, idx].imshow(stack[:,:,end - n] .- ref, vmax=dvm, vmin=-dvm, cmap=cmap)
    end

    return fig, axes
end
