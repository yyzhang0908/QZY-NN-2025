function pe_sig = get_peri_event_signal(sig, ts, ev, twin)
itpt = ev(:) + twin(:)';
pe_sig = interp1(ts, sig, itpt);
end